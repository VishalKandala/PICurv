# Smoke Fixtures Guide

This directory holds lightweight executable smoke assets for the canonical `make smoke` target.

The current smoke runner focuses on entrypoint-level binary validation:

- `bin/simulator` must launch and respond to `-help`
- `bin/postprocessor` must launch and respond to `-help`

The subdirectories below are reserved for future tiny runtime fixtures that drive true case-level
smoke runs without relying on the large example suite.
