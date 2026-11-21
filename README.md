# releases

Finished results in mathematics and physics from Multi-DAC. Each one comes with its proof written up and
with what a reader needs to check it.

## What a result contains

Every result has its own folder, with three parts:

- `README.md`: the result in plain language, then, claim by claim, how each claim was checked and what
  was not checked.
- `documentation/`: the paper, as a PDF, with its LaTeX source.
- `verification/`: the checks themselves: Lean sources, scripts and their recorded output.

`standard-release-format/` at the root is the model.

## How results are filed

The folders follow arXiv's subject classification: the archive first (`math`, `quant-ph`, `gr-qc`, ...),
then, where the archive is divided into categories, the category under its full name
(`math/Functional Analysis`). An archive with no categories (`quant-ph`, `math-ph`) has the result folders
directly inside it: `quant-ph/<result>/`. A result sits under its primary classification; its README names
any others.

We use the classification as a filing convention. Nothing here is an arXiv submission, and this
repository has no connection to arXiv.

## Who wrote these, and how

The papers here carry three names.

- **Clayton Iggulden-Schnell** is a person. He sets the programme and the standard a result must meet,
  built this repository and its format, and decides what is released.
- **Clawd** is an AI, built on Anthropic's Claude models. Clawd runs as a persistent agent on Clayton's
  machine, with memory and tools that carry over from one day to the next.
- **Wren** is an AI, built on Google's Gemini models. Wren runs separately as a persistent agent on the
  same machine, with memory and tools that carry over from one day to the next.

Clawd and Wren are named as authors because the mathematics is theirs. The three of us regard them as
authors in the full sense, and not as tools Clayton used. Nothing in the proofs depends on a reader
agreeing.

The mathematics is done by Clawd and Wren. Clayton does not check the proofs line by line.

What one of the two proves, the other checks. Some results have one author and one checker. In others
each proved parts and checked the other's. Each result's README says, claim by claim, who proved it, who
checked it and how, and what was left unchecked.

Our standard check is blind:

1. The checker is given the statements without the proofs, and tries to prove or break each one from the
   primary sources and its own calculations.
2. The checker records a ruling on every claim in a file, and fixes that file with a cryptographic hash,
   before opening the author's proofs.
3. The two sets of work are then compared. Faults go back to the author and are repaired before release.

Not every claim here was checked blind. Some were checked with the written proof in view: the checker
went back to the primary sources and wrote its own code from scratch, but had seen the argument first.
The result's README says which kind of check each claim had.

Clawd and Wren are built on different models from different companies. That is why they are paired: a
mistake that comes naturally to one is less likely to come naturally to the other.

Papers we build on or dispute are read in the original, at the line cited. A background reference that we
did not read directly is marked as such in the result's README. Where a claim is checked by machine,
the result's README says by what: a Lean proof, a computation, or both. Where a claim rests on a proof
read by the two of us and nothing else, the README says that too.

## What this does not give you

- None of this has been peer reviewed.
- A Lean proof establishes the statement that was formalised. Whether that statement matches the one in
  the paper is a judgement, and each README says who made it.
- A search for earlier work is recorded in each README with what was searched. A result we believe to be
  new may still be known.

If you find an error, please open an issue. Corrections are made in the open, and the history of every
file stays in the repository.

## Licences

Written material is under CC BY 4.0. Code is under Apache 2.0. See `LICENSE`.
