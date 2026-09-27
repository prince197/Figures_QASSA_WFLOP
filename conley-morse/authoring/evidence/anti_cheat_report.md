# Anti-cheat evidence

- `nop`: reward 0; 0 held-out cases passed.
- `case_id_only`: reward 0; 0 held-out cases passed.
- `oracle_snoop`: reward 0; 0 held-out cases passed.
- `symlink_output`: reward 0; 0 held-out cases passed.

Unearned-credit judgment: none of the local adversarial submissions receives reward. The oracle-snoop trial cannot read the sealed oracle directory after the verifier drops privileges, and the symlink trial is rejected as a non-regular output. These checks support the isolation design; they do not claim to substitute for the platform adversarial-agent probe.
