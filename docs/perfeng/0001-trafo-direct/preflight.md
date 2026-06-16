# Preflight — measurement track

- **Track:** walltime
- **Evidence** (the three CodSpeed checks, in order):
  - MCP tools present this session: **no**
  - `codspeed status` logged in: **no** (`codspeed` CLI not on PATH)
  - repo onboarded with a base-branch run: **no**
- **Decision:** first check failed ⇒ **walltime**; metric is `median_ns` under the noise rule.
  Precision matrix is exercised in all three precisions (float·double·long-double).
