# Security Policy

## Scope

Blind Signal Discovery (BSD) is a standalone local research prototype intended for offline signal-processing and statistical experiments on synthetic and locally supplied time-series data. It is not designed or hardened as a public multi-tenant web service.

- The optional Streamlit Workstation (`app/bsd_workstation.py`) binds to `127.0.0.1` by default and should only be run locally on trusted files.
- Do not expose the local Workstation port (`8501`) directly to untrusted networks.
- Third-party public datasets (such as GWOSC / LIGO Open Science Center strain releases) remain subject to their respective upstream data licenses and terms of use.

## Reporting a Vulnerability

If you discover a security issue in this repository:

1. Please report it privately via **GitHub Security Advisories** ("Report a vulnerability" under the **Security** tab of this repository), or open a minimal issue without publishing exploit details if private reporting is unavailable.
2. Include a concise description of the issue, reproduction steps, and affected files or commit hashes.
3. Please allow reasonable time for triage and remediation before public disclosure.
