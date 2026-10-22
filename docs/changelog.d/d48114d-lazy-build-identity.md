- Faster CLI on network filesystems: `picurv --help`, `validate` and `init` no longer run git, the
  identity takes one `git describe` instead of three or four git calls, and `picurv version` reads each
  executable's identity from a stamp in its file instead of starting it. Rebuild to add the stamp.
