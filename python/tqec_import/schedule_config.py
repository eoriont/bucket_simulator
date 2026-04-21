from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class RoundSchedule:
    pre_rounds: int
    merge_rounds: int
    post_rounds: int


def circuit_distance_from_k(circuit_k: int) -> int:
    return 2 * circuit_k + 1


def default_round_schedule_for_k(circuit_k: int) -> RoundSchedule:
    d = circuit_distance_from_k(circuit_k)
    return RoundSchedule(
        pre_rounds=d + 1,
        merge_rounds=d + 1,
        post_rounds=0,
    )


def round_schedule_to_dict(schedule: RoundSchedule) -> dict[str, int]:
    return asdict(schedule)
