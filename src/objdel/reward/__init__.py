"""Rewards built on the checker: scoring against an annotator set and GRPO group advantages."""
from .annotators import AnnotatorSet, Score, group_advantages

__all__ = ["AnnotatorSet", "Score", "group_advantages"]
