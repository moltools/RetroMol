from pathlib import Path

import numpy as np


TOKEN_UNK = "<UNK>"
SPECIAL_TOKENS = [TOKEN_UNK]


class Vocabulary:
    def __init__(self, tokens: list[str] | set[str]) -> None:
        ordered: list[str] = []

        for token in SPECIAL_TOKENS:
            ordered.append(token)

        for token in sorted(tokens):
            if token not in ordered:
                ordered.append(token)

        self.tokens = ordered
        self.stoi = {token: i for i, token in enumerate(self.tokens)}
        self.itos = {i: token for i, token in enumerate(self.tokens)}

    def __len__(self) -> int:
        return len(self.tokens)

    @property
    def unk_id(self) -> int:
        return self.stoi[TOKEN_UNK]

    @classmethod
    def from_tokens(cls, tokens: list[str]) -> "Vocabulary":
        return cls(tokens)

    @classmethod
    def from_file(cls, path: str | Path) -> "Vocabulary":
        path: Path = Path(path)
        tokens = path.read_text().splitlines()
        return cls.from_tokens(tokens)

    def write(self, path: str | Path) -> None:
        path: Path = Path(path)
        path.write_text("\n".join(self.tokens))


class Fingerprinter:
    def __init__(self, vocab: Vocabulary, n_bits: int | None = None, n_hashes: int = 1) -> None:
        if n_bits is not None and n_bits < 1:
            raise ValueError("n_bits must be >= 1")

        if n_hashes < 1:
            raise ValueError("n_hashes must be >= 1")

        self.vocab = vocab
        self.n_bits = n_bits
        self.n_hashes = n_hashes

    @property
    def size(self) -> int:
        if self.n_bits is None:
            return len(self.vocab)

        return self.n_bits

    @property
    def should_fold(self) -> bool:
        if self.n_bits is None:
            return False

        return len(self.vocab) > self.n_bits

    def _indices(self, token_id: int) -> list[int]:
        if not self.should_fold:
            return [token_id]

        indices = []

        for i in range(self.n_hashes):
            # Deterministic integer hash
            x = token_id + i * 0x9E3779B9
            x = (x ^ (x >> 16)) * 0x85EBCA6B
            x = (x ^ (x >> 13)) * 0xC2B2AE35
            x = x ^ (x >> 16)

            indices.append(x % self.size)

        return list(dict.fromkeys(indices))

    def encode(self, per_monomer_tokens: list[list[str]]) -> np.ndarray:
        # Every sublist represents tokens representing a single monomer (i.e., name + higher level pseudonyms)
        fp = np.zeros((self.size,), dtype=np.float32)

        for tokens in per_monomer_tokens:
            # Avoid duplicate tokens within same monomer
            tokens = list(dict.fromkeys(tokens))

            if not tokens:
                indices = self._indices(self.vocab.unk_id)
                weight = np.float32(1.0 / len(indices))

                for index in indices:
                    fp[index] += weight

                continue

            weight = np.float32(1.0 / len(tokens))

            for token in tokens:
                token_id = self.vocab.stoi.get(token, self.vocab.unk_id)
                indices = self._indices(token_id)
                per_bin_weight = np.float32(weight / len(indices))

                for index in indices:
                    fp[index] += per_bin_weight

        return fp

    @staticmethod
    def similarity(fp1: np.ndarray, fp2: np.ndarray) -> float:
        # Cosine similarity
        dot_product = float(np.dot(fp1, fp2))
        norm_fp1 = float(np.linalg.norm(fp1))
        norm_fp2 = float(np.linalg.norm(fp2))

        if norm_fp1 == 0.0 or norm_fp2 == 0.0:
            return 0.0

        return dot_product / (norm_fp1 * norm_fp2)
