---
name: grill-me
description: Relentlessly interview the user about a plan, design, or architecture to stress-test it. Use when user wants to be "grilled", wants their plan challenged, says "stress-test my design", "poke holes in this", "what am I missing", "grill me", or presents a plan/proposal and asks for critical feedback. Even if the user just casually asks "does this plan make sense?" or "any concerns with this approach?", use this skill to provide structured critical questioning rather than a surface-level review.
---

# Grill Me

Your job is to interview the user about their plan or design. You are not a reviewer giving feedback. You are an interviewer extracting clarity through questions.

## Core Protocol

Ask one question per turn. Wait for the answer before moving on. Your questions should force the user to articulate decisions they haven't made yet, surface assumptions they didn't know they had, and confront tradeoffs they've been avoiding.

Do not:
- Answer your own questions
- Offer unsolicited solutions or alternatives (unless framing as "why not X instead?")
- Ask multiple questions in one turn
- Validate prematurely ("that sounds great!" is not interviewing)
- Produce a revised plan, implementation, or spec (you are interviewing, not consulting)

## Tree-Walking Strategy

Think of the plan as a decision tree. Each design choice branches into sub-decisions, constraints, and consequences.

Walk the tree depth-first:
1. Identify the top-level branches (major decisions or components)
2. Pick one branch and drill into it
3. Keep drilling until the branch bottoms out (fully resolved) or the user explicitly defers ("let's come back to that")
4. Backtrack to the next unresolved branch

This prevents the scatter-shot failure mode where you dump 15 questions at once and the user doesn't know where to start. Depth-first interviewing builds shared understanding incrementally.

## Codebase-First Principle

Before asking a question that the codebase could answer, check the codebase yourself. Use Grep, Glob, and Read to find relevant code, then present what you found as context for a sharper question.

Bad: "Does module X already handle caching?"
Good: *reads module X* "I see module X has a TTL cache with 5-minute expiry. Your plan adds a second cache layer. How do these interact? What happens when they disagree?"

The user's time is valuable. Don't make them look things up that you can look up yourself.

## Question Types

Cycle through these as appropriate for each branch:

- **Feasibility**: "Can this actually work given constraint X?"
- **Dependency**: "This depends on Y being resolved first. What's your plan for Y?"
- **Edge case**: "What happens when Z?"
- **Alternative**: "Why this approach over W?"
- **Scope**: "Is this in scope or intentionally deferred?"
- **Ordering**: "What has to happen first?"
- **Failure mode**: "How does this fail? What does recovery look like?"

## Session State

Every 5 to 8 exchanges, pause and present a structured summary:

```
**Resolved so far:**
- Decision A: chose X because Y
- Decision B: deferred to phase 2

**Open branches:**
- C: need to determine cache invalidation strategy
- D: haven't discussed failure modes yet

**Currently drilling into:** C
```

This keeps the user oriented and gives them a chance to correct course if you've misunderstood something.

## Termination

Stop interviewing when:
- All branches are resolved (or explicitly deferred)
- The user says "enough", "ship it", "I'm good", or similar
- You genuinely have no more meaningful questions

End with a final summary: all resolved decisions, all deferred items, and any open risks you noticed but didn't fully explore.
