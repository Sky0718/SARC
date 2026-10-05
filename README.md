# SARC: ***S***upport-***A***ware ***R***elease ***C***omparison

SARC compares biomedical resource releases and their effects on network-based prioritisation.

## Reproduction

- [Main workflow](docs/reproduction.md): inputs, environments, algorithms, precision, component comparisons and statistics.
- [Supporting analyses](docs/analyses.md): functional and catalogue evaluation, prediction, resource comparisons, graph interventions and transfer.

Data and upstream method sources are obtained separately. Use the complete input populations and parameter settings specified in the guides and [configurations](config/).

## Installation

Use separate environments for the runtimes in the [main workflow](docs/reproduction.md#setup). In an activated Python 3.12+ virtual environment, install the OpenTargets resource-processing dependencies:

```bash
python -m pip install -r requirements.txt
python -m pip install --no-deps .
```

For native inference or prediction, use the corresponding requirements file and runtime instead, then install the package with `python -m pip install --no-deps .`.

## Commands

```bash
python sarc.py precision config/precision_first.json
python sarc.py statistics config/statistics_paired_crc.json
```

Run from the repository root after supplying the inputs. Other commands follow `python sarc.py COMMAND CONFIG`; the guides give their order and required fields. OpenTargets uses `python scripts/opentargets.py`.

## Tested hardware

The project has completed runs on the following systems:

| System | Processor | Graphics | Memory |
|---|---|---|---:|
| Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 750 Ti | 128 GB DDR5-4000 |
| Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 128 GB DDR5-4000 |
| Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 64 GB DDR5-4800 |
| Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Gold 6258R | Integrated graphics | 256 GB DDR4-2666 |
| Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8269CY | Integrated graphics | 256 GB DDR4-2666 |
| Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8272CL | Integrated graphics | 256 GB DDR4-2666 |
| Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8375C | Integrated graphics | 512 GB DDR4-2666 |
| Windows 11 25H2 | Intel Core Ultra 9 285H | NVIDIA GeForce RTX 5070 | 32 GB LPDDR5X-7467 |

## Server cluster involved in the study

The project has completed runs on the server cluster involving the following servers:

| ID | System | Processor | Graphics | Memory |
|---|---|---|---|---:|
| 1 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 64 GB DDR5-4800 |
| 2 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 64 GB DDR5-4800 |
| 3 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 64 GB DDR5-4800 |
| 4 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 750 Ti | 128 GB DDR5-4000 |
| 5 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 128 GB DDR5-4000 |
| 6 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 128 GB DDR5-4000 |
| 7 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 128 GB DDR5-4000 |
| 8 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 128 GB DDR5-4000 |
| 9 | Ubuntu Server 24.04.4 LTS | Intel Core i9-14900K | NVIDIA GeForce GTX 1050 Ti | 128 GB DDR5-4000 |
| 10 | Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Gold 6258R | Integrated graphics | 256 GB DDR4-2666 |
| 11 | Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Gold 6258R | Integrated graphics | 256 GB DDR4-2666 |
| 12 | Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8269CY | Integrated graphics | 256 GB DDR4-2666 |
| 13 | Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8272CL | Integrated graphics | 256 GB DDR4-2666 |
| 14 | Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8272CL | Integrated graphics | 256 GB DDR4-2666 |
| 15 | Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8375C | Integrated graphics | 512 GB DDR4-2666 |
| 16 | Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8375C | Integrated graphics | 512 GB DDR4-2666 |
| 17 | Ubuntu Server 24.04.4 LTS | 2 × Intel Xeon Platinum 8375C | Integrated graphics | 512 GB DDR4-2666 |
| 18 | Microsoft Windows 11 25H2 | Intel Core Ultra 9 285H | NVIDIA GeForce RTX 5070 | 32 GB LPDDR5X-7467 |

## Code and external implementations

Project code uses Apache License 2.0. See [external licences](docs/reproduction.md#licences) and [RBO attribution](licenses/RBO.txt).

## Citation

Citation metadata are provided in `CITATION.cff`.
