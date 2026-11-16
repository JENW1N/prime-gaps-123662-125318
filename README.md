# Consecutive primes at distances 123,662 and 125,318

**Will Howe, Grant Huffman, Jensen Kohlmeyer, and Liam Kruer**

[Download the proof paper](https://github.com/JENW1N/prime-gaps-123662-125318/releases/latest/download/consecutive-primes-123662-125318-proof.pdf)

The four-page paper gives computer-assisted proofs of two exact consecutive-prime gaps. It contains two pages of mathematical argument, two pages displaying the starting primes, and the complete arithmetic appendix embedded in the PDF. OpenAI Codex is acknowledged for assistance with construction, arithmetic verification, and preparation.

| Gap | Endpoint digits | Covering moduli | Largest covering prime | ECPP steps, left and right |
| --- | ---: | ---: | ---: | --- |
| 123,662 | 8,416 | 2,216 | 32,587 | 679 and 671 |
| 125,318 | 6,454 | 1,754 | 14,983 | 538 and 537 |

For 123,662, covering congruences certify all 123,661 interior integers composite. For 125,318, the congruences cover 124,526 interior integers, and the remaining 791 have 175 proper-divisor witnesses and 616 base-2 Fermat compositeness witnesses.

Both full verification replays completed successfully on November 16, 2026. They prove all four endpoints prime and check every interior integer. The replay transcripts are included.

## Verification

With Python and gmpy2 installed, run these commands from the repository directory.

```sh
python verify_gap_122290.py --gap 123662 --left gap_123662_left.ecpp --right gap_123662_right.ecpp --jobs 8
python verify_gap_122290.py --gap 125318 --left gap_125318_left.ecpp --right gap_125318_right.ecpp --jobs 8
```

The verifier reconstructs the endpoints from the CRT data, checks the interior witnesses, and verifies the elliptic-curve certificate chains. The paper also supplies elementary order certificates for the four terminal primes.

## Publication and record status

The paper and complete certificates are public here. The source archive in the release supplies the mathematical appendix in twelve parts for Hexagon's per-file size limit.

Prime Gap Records submissions were attempted on November 16, 2026. Its server deferred them because of its recent-submission limit, following the accepted 122,290 record. [Catalog submission issue #40](https://github.com/primegap-list-project/prime-gap-list/issues/40) includes both records and asks the maintainers to review or import them. Hexagon received the joint paper as **hexagon:2611.00536** on November 16, 2026. Its preview compiled successfully. Processing and moderator approval are pending. [Submission receipt](https://hexagonmath.org/submit/7e4e15b1-a98f-4b2c-a6ed-2c8d54d3f510) requires the submitting account.

Both gap sizes are absent from the catalog snapshot of November 29, 2019 and from the catalog checked on November 16, 2026. Their original problem dates have not been established.

[Historical catalog](https://github.com/primegap-list-project/primegap-list-project.github.io/blob/4a12d01a11bd48daa37a0fedc4e7e0cf564d950e/_data/allgaps.csv)

[Earlier certified gap of 122,290](https://github.com/JENW1N/prime-gap-122290)
