// Actual Chromium checks: file-origin loading, non-decision safety, accessible filtering and independent scroll.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const {pathToFileURL} = require('node:url');
const {chromium} = require('playwright');

(async () => {
  const root = path.resolve(process.argv[2]);
  const browser = await chromium.launch({headless: true, args: ['--no-sandbox']});
  const results = [];
  try {
    for (const mode of ['compact_single_file', 'paged_single_file', 'data_backed']) {
      for (const state of ['resolved', 'unresolved']) {
        const page = await browser.newPage({viewport: {width: 1400, height: 1000}});
        const errors = [], external = [];
        page.on('pageerror', e => errors.push(e.message));
        page.on('request', r => {if (/^https?:/.test(r.url())) external.push(r.url());});
        const file = path.join(root, `${mode}-${state}.html`);
        await page.goto(pathToFileURL(file).href);
        if (mode === 'compact_single_file') {
          await page.waitForSelector('#idx tbody tr');
          assert.equal(await page.locator('#evidence-unresolved').count(), state === 'resolved' ? 0 : 1);
          assert.equal(await page.locator('#n_clean_pass').innerText(), state === 'resolved' ? '300' : '150');
          const before = await page.evaluate(() => JSON.stringify({gates:G.map(g=>g.gates), thresholds:T}));
          if (state === 'unresolved') {
            await page.locator('#evidence-unresolved > summary').click();
            await page.getByRole('button', {name: 'Show unresolved guides'}).focus();
            await page.keyboard.press('Enter');
            assert.equal(await page.locator('#idx tbody tr').count(), 150);
            assert.equal(await page.locator('[data-preset="unresolved"]').getAttribute('aria-pressed'), 'true');
            await page.locator('#idx tbody tr').first().click();
            assert.equal(await page.locator('#detail [data-warning-id="GATE_RUN_UNAVAILABLE"]').count(), 1);
            assert.equal(await page.evaluate(() => JSON.stringify({gates:G.map(g=>g.gates), thresholds:T})), before);
            await page.reload();await page.waitForSelector('#idx tbody tr');
            assert.equal(await page.locator('#idx tbody tr').count(),150);
            assert.equal(await page.locator('[data-preset="unresolved"]').getAttribute('aria-pressed'),'true');
          }
          if (!(await page.locator('#filters').evaluate(e=>e.open))) await page.locator('#filters > summary').click();
          await page.waitForTimeout(100);
          const scroll = await page.evaluate(() => {
            const filters=document.querySelector('#filters'), index=document.querySelector('#index-scroll');
            const top=filters.getBoundingClientRect().top;
            index.scrollTop=1000;
            const filterBefore=filters.scrollTop, indexBefore=index.scrollTop;
            filters.scrollTop=1000;
            return {filterOverflow:filters.scrollHeight>filters.clientHeight,
              indexOverflow:index.scrollHeight>index.clientHeight,
              indexBefore,indexAfter:index.scrollTop,filterBefore,filterAfter:filters.scrollTop,
              fixed:Math.abs(filters.getBoundingClientRect().top-top)<1};
          });
          assert(scroll.filterOverflow && scroll.indexOverflow, JSON.stringify(scroll));
          assert(scroll.indexBefore > 0 && scroll.filterAfter > scroll.filterBefore && scroll.fixed, JSON.stringify(scroll));
          assert.equal(scroll.indexAfter, scroll.indexBefore);
          const headerOffset=await page.evaluate(()=>document.querySelector('#idx th').getBoundingClientRect().top-document.querySelector('#index-scroll').getBoundingClientRect().top);
          assert(Math.abs(headerOffset)<2, `table header displaced within index: ${headerOffset}`);
          const quiet = page.locator('#frows > details > summary');
          if (await quiet.count()) await quiet.click();
          await page.locator('.gate-meta > summary').first().scrollIntoViewIfNeeded();
          await page.locator('.gate-meta > summary').first().click();
          const metadata = await page.locator('.gate-meta').first().innerText();
          for (const label of ['Setting key', 'Run comparison', 'Action', 'Stage', 'Scope', 'Units', 'Definition']) assert(metadata.includes(label), label);
          await page.setViewportSize({width:800,height:900});
          assert(await page.evaluate(() => document.querySelector('#index-scroll').clientHeight>0));
        } else {
          await page.waitForSelector('#rows tr');
          assert.equal(await page.locator('#evidence-notices').count(), state === 'resolved' ? 0 : 1);
          const before = await page.evaluate(() => JSON.stringify(G.map(g=>g.gates)));
          if (state === 'unresolved') {
            await page.getByLabel('Show unresolved evidence (150)').check();
            assert((await page.locator('#visible-count').innerText()).startsWith('150 guides'));
            await page.locator('#rows button[data-guide]').first().click();
            await page.waitForSelector('#details [data-warning-id="GATE_RUN_UNAVAILABLE"]');
            assert.equal(await page.evaluate(() => JSON.stringify(G.map(g=>g.gates))), before);
          }
          assert.equal(await page.getByRole('complementary', {name:'Evidence and thresholds'}).count(), 1);
          await page.locator('aside details:has(> dl) > summary').first().click();
          const scroll=await page.evaluate(()=>{
            const aside=document.querySelector('main>aside'),index=document.querySelector('main>section');
            const top=aside.getBoundingClientRect().top;
            index.scrollTop=800;const indexBefore=index.scrollTop;
            aside.scrollTop=800;
            return {indexBefore,indexAfter:index.scrollTop,aside:aside.scrollTop,
              fixed:Math.abs(aside.getBoundingClientRect().top-top)<1};
          });
          assert(scroll.indexBefore>0&&scroll.aside>0&&scroll.fixed,JSON.stringify(scroll));
          assert.equal(scroll.indexBefore,scroll.indexAfter);

        }
        const overflow=await page.evaluate(()=>document.body.scrollHeight-innerHeight);
        assert(overflow<=1, `page exceeds viewport by ${overflow}px: ${mode}/${state}`);
        assert.deepEqual(errors, []);
        assert.deepEqual(external, []);
        await page.screenshot({path:path.join(root, `${mode}-${state}.png`),fullPage:true});
        results.push({mode,state,passed:true});
        await page.close();
      }
    }
    fs.writeFileSync(path.join(root, 'browser-results.json'), JSON.stringify(results, null, 2)+'\n');
    console.log(JSON.stringify(results));
  } finally { await browser.close(); }
})().catch(e => {console.error(e); process.exitCode=1;});
