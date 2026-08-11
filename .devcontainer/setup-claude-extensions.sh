#!/usr/bin/env bash
#
# Install a well-defined set of Claude Code plugins and skills into the dev
# container, and apply user-level Claude Code settings, so a freshly started
# container has them by default.
#
# Scope is user-level (~/.claude). The script is idempotent and never aborts
# container creation on a failed install.
#
# What it installs:
#   Plugins:
#     - superpowers@claude-plugins-official   (Anthropic's Superpowers)
#     - codspeed@claude-plugins-official       (CodSpeed)
#   Skills:
#     - mattpocock/skills                      (Matt Pocock's skills)
#     - shadcn/improve                         (shadcn Improve skill)
#   Settings (~/.claude/settings.json):
#     - attribution.commit = ""                (no Claude attribution on commits)
#     - attribution.pr     = ""                (no Claude attribution on PRs)

set -u

MARKETPLACE="claude-plugins-official"
PLUGINS=(
  "superpowers@${MARKETPLACE}"
  "mattpocock-skills@${MARKETPLACE}"
)
# GitHub owner/repo[/subpath] shorthand resolved by the `skills` CLI.
SKILL_SOURCES=()

failures=0

log()  { printf '[claude-extensions] %s\n' "$*"; }
warn() { printf '[claude-extensions] WARNING: %s\n' "$*" >&2; failures=$((failures + 1)); }

# Run a command; log ✓/✗ but never abort the script.
try() {
  local label="$1"; shift
  if "$@"; then
    log "✓ ${label}"
  else
    warn "✗ ${label} (exit $?) — continuing; re-run this script to retry"
  fi
}

# --- Preflight -------------------------------------------------------------
# Tooling is provided by the claude-code devcontainer feature. If it's absent,
# the feature ordering is wrong; warn loudly but don't fail container creation.
if ! command -v claude >/dev/null 2>&1; then
  warn "'claude' CLI not found on PATH — skipping plugin install."
  CLAUDE_OK=0
else
  CLAUDE_OK=1
fi
if ! command -v npx >/dev/null 2>&1; then
  warn "'npx' not found on PATH — skipping skills install."
  NPX_OK=0
else
  NPX_OK=1
fi

# --- Plugins ---------------------------------------------------------------
if [ "${CLAUDE_OK}" = "1" ]; then
  # The official marketplace is auto-known, but ensure it's registered so a
  # pristine ~/.claude still resolves plugin@marketplace. `marketplace add` on an
  # already-known marketplace is a no-op-ish error, so tolerate it.
  if ! claude plugin marketplace list 2>/dev/null | grep -q "${MARKETPLACE}"; then
    try "marketplace add ${MARKETPLACE}" \
      claude plugin marketplace add "anthropics/${MARKETPLACE}"
  fi

  # Snapshot already-installed plugins to keep re-runs quiet.
  installed_plugins="$(claude plugin list 2>/dev/null || true)"
  for plugin in "${PLUGINS[@]}"; do
    if printf '%s\n' "${installed_plugins}" | grep -q "${plugin%@*}"; then
      log "• ${plugin} already installed — skipping"
      continue
    fi
    try "install ${plugin}" \
      claude plugin install "${plugin}" --scope user
  done
fi

# --- Skills ----------------------------------------------------------------
if [ "${NPX_OK}" = "1" ]; then
  for src in "${SKILL_SOURCES[@]}"; do
    # --skill '*'        : take every skill in the source repo
    # --agent claude-code: deterministic target (don't rely on auto-detection)
    # --global           : install into ~/.claude/skills
    # --yes              : accept defaults, no prompts
    try "skills add ${src}" \
      npx -y skills@latest add "${src}" \
        --skill '*' --agent claude-code --global --yes
  done
fi

# --- Settings --------------------------------------------------------------
# Apply user-level Claude Code settings. Merge into any existing settings.json
# rather than clobbering it, preserving keys the user (or other tooling) set.
# `node` ships with the claude-code devcontainer feature; fall back gracefully.
if command -v node >/dev/null 2>&1; then
  # Empty attribution strings = no "Co-Authored-By: Claude" / "Generated with
  # Claude Code" text on commits or PRs. (Supersedes the deprecated
  # `includeCoAuthoredBy` boolean.)
  try "clear git attribution" node -e '
    const fs = require("fs");
    const os = require("os");
    const path = require("path");
    const dir = path.join(os.homedir(), ".claude");
    const file = path.join(dir, "settings.json");
    fs.mkdirSync(dir, { recursive: true });
    let settings = {};
    try {
      settings = JSON.parse(fs.readFileSync(file, "utf8"));
    } catch (e) { /* missing or empty/invalid: start fresh */ }
    settings.attribution = { ...settings.attribution, commit: "", pr: "" };
    fs.writeFileSync(file, JSON.stringify(settings, null, 2) + "\n");
  '
else
  warn "'node' not found on PATH — skipping settings.json update."
fi

# --- Summary ---------------------------------------------------------------
if [ "${failures}" -eq 0 ]; then
  log "All plugins and skills installed."
else
  log "Finished with ${failures} issue(s); container creation not blocked. Re-run to retry."
fi

# Always succeed: a transient install failure must not break the dev environment.
exit 0
