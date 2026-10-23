# Temporal graph benchmark integration

The canonical E47-style projection research and reported synthetic link-prediction experiments live in [E47-Kartekeya / temporal-graph-learning](https://github.com/nicholaskouns-create/E47-Kartekeya/tree/main/temporal-graph-learning).

This repository remains the canonical implementation home for the quantum clock, proper-time dynamics and TQGNN engine. The benchmark scaffold uses a random rank-k orthogonal projector matched to occupancy 47/125; it is **not** a direct implementation of the exact E47 Casimir projector. No published A2QTGN or official TGBL result has been reproduced by that scaffold.

**Integration boundary:** expose node embeddings from this repo's temporal model, insert the projection on the *feature axis*, then apply the existing link predictor. Validate on chronological splits with identical negative sampling and multiple seeds. See the canonical benchmark README for reported metrics and caveats.
