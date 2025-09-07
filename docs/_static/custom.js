// Custom JavaScript for siRNAforge documentation

// Load Mermaid library and initialize diagrams
document.addEventListener('DOMContentLoaded', function() {
    // Load Mermaid from CDN
    var script = document.createElement('script');
    script.src = 'https://cdn.jsdelivr.net/npm/mermaid@10.6.1/dist/mermaid.min.js';
    script.onload = function() {
        // Initialize Mermaid with custom theme
        mermaid.initialize({
            startOnLoad: true,
            theme: 'default',
            themeVariables: {
                primaryColor: '#2980b9',
                primaryTextColor: '#2c3e50',
                primaryBorderColor: '#3498db',
                lineColor: '#34495e',
                secondaryColor: '#ecf0f1',
                tertiaryColor: '#e74c3c',
                background: '#ffffff',
                mainBkg: '#ffffff',
                secondBkg: '#f8f9fa'
            },
            flowchart: {
                useMaxWidth: true,
                htmlLabels: true,
                curve: 'basis'
            },
            sequence: {
                useMaxWidth: true,
                wrap: true,
                diagramMarginX: 50,
                diagramMarginY: 10
            },
            gantt: {
                useMaxWidth: true
            }
        });
        
        // Process all mermaid diagrams
        mermaid.run();
    };
    document.head.appendChild(script);
    
    // Enhanced search functionality
    function enhanceSearch() {
        var searchInput = document.querySelector('input[type="search"]');
        if (searchInput) {
            searchInput.setAttribute('placeholder', '🔍 Search siRNAforge docs...');
        }
    }
    
    // Call enhancements
    enhanceSearch();
    
    // Add scroll-to-top functionality
    var scrollButton = document.createElement('button');
    scrollButton.innerHTML = '↑';
    scrollButton.className = 'scroll-to-top';
    scrollButton.style.cssText = `
        position: fixed;
        bottom: 20px;
        right: 20px;
        background: #2980b9;
        color: white;
        border: none;
        border-radius: 50%;
        width: 50px;
        height: 50px;
        cursor: pointer;
        display: none;
        z-index: 1000;
        font-size: 18px;
        box-shadow: 0 2px 5px rgba(0,0,0,0.2);
    `;
    
    scrollButton.addEventListener('click', function() {
        window.scrollTo({top: 0, behavior: 'smooth'});
    });
    
    document.body.appendChild(scrollButton);
    
    // Show/hide scroll button
    window.addEventListener('scroll', function() {
        if (window.pageYOffset > 300) {
            scrollButton.style.display = 'block';
        } else {
            scrollButton.style.display = 'none';
        }
    });
});
