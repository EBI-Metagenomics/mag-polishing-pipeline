Empty placeholder files.

`conf/test.config` needs paths that satisfy `exists: true` in `assets/schema_input.json`,
and the `test` profile is a DAG-resolution smoke test run with `-preview` or `-stub-run`,
so nothing here is ever read. A real test profile needs a tiny eukaryotic metagenome plus a
Branchwater index, neither of which is small enough to ship.
