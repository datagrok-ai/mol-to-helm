After every code change:
1. Run all tests (Linear, Cyclic, BILN). Baselines: 20/20, 41/41, 1/5. No regressions allowed.
2. Update docs (README.md, DESCRIPTION.md, PROJECT_CONTEXT.md, CLAUDE.md) if behavior changed.
3. Run `cd aggregated && py aggregate_logics.py` if any logics/ file changed. Verify syntax of output.
