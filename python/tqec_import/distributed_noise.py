from __future__ import annotations

from dataclasses import asdict, dataclass
from math import exp, isnan, pow
from pathlib import Path
from typing import Iterable

import stim


DISTILLATION_PROTOCOLS = (
    "none",
    "radar",
    "pumping_2to1",
    "pumping_3to1",
    "recurrence_2to1",
    "recurrence_3to1",
)
REMOTE_2Q_GATES = ("CX", "CNOT", "CZ")
LOCAL_2Q_GATES = ("CX", "CNOT", "CZ")  # all 2Q gates get physical depolarizing noise


@dataclass
class NoiseConfig:
    physical_error: float = 0.001
    interconnect_error: float = 0.0
    accurate_rcx: bool = False
    channel_depolarization_error: float | None = None
    raw_epr_fidelity: float = 0.99
    distillation_protocol: str = "none"
    distillation_rounds: int = 1
    distillation_backup_batches: int = 1
    entanglement_rate: float = 100e6
    measurement_delay: float = 660e-9
    t1_coherence_time: float = 125e-6
    t2_coherence_time: float = 200e-6
    monolithic_baseline: bool = False
    local_baseline: bool = False
    split_axis: str = "auto"
    split_x: float | None = None


@dataclass
class DistillationResult:
    protocol: str
    output_fidelity: float
    success_probability: float
    raw_pairs_consumed: int
    distillation_time: float


@dataclass
class SectionNoiseSummary:
    section_index: int
    repeat_count: int
    remote_cnots_per_body: int
    required_epr_pairs_per_round: int
    effective_distillation_protocol: str
    raw_channel_error: float
    effective_channel_error: float
    distilled_fidelity: float
    distilled_error: float
    distillation_success_probability: float
    probability_all_distillation_fail: float
    remote_cnot_error: float
    distillation_time_ns: float
    idling_time_us: float
    p_x: float
    p_y: float
    p_z: float
    timing_constraint_satisfied: bool
    raw_pairs_per_distilled: int
    distillation_backup_batches: int
    distillation_epr_slots_required: int
    distillation_qubits_required: int


def protocol_branching_factor(protocol: str) -> int:
    if protocol in ("pumping_3to1", "recurrence_3to1"):
        return 3
    if protocol in ("pumping_2to1", "recurrence_2to1"):
        return 2
    return 1


def is_recurrence_protocol(protocol: str) -> bool:
    return protocol in ("recurrence_2to1", "recurrence_3to1")


def is_pumping_protocol(protocol: str) -> bool:
    return protocol in ("pumping_2to1", "pumping_3to1")


def compute_2to1_success_prob(f1: float, f2: float) -> float:
    term1 = (1.0 - 2.0 * (1.0 - f1) / 3.0) * (1.0 - 2.0 * (1.0 - f2) / 3.0)
    term2 = 4.0 * (1.0 - f1) * (1.0 - f2) / 9.0
    return max(0.01, term1 + term2)


def compute_2to1_fidelity(f1: float, f2: float) -> float:
    p_s = compute_2to1_success_prob(f1, f2)
    if p_s <= 0:
        return f2
    numerator = (f1 * f2) + ((1.0 - f1) * (1.0 - f2) / 9.0)
    return min(1.0, max(0.0, numerator / p_s))


def compute_3to1_success_prob(f1: float, f2: float) -> float:
    f1_sq = f1 * f1
    one_minus_f1 = 1.0 - f1
    p_s = (
        f1_sq * f2
        + f1_sq * (1.0 - f2) / 9.0
        + 2.0 * f1 * one_minus_f1 * (f2 / 3.0 + (1.0 - f2) / 9.0)
        + one_minus_f1 * one_minus_f1 * (f2 / 9.0 + (1.0 - f2) / 81.0)
    )
    return max(0.01, p_s)


def compute_3to1_fidelity(f1: float, f2: float) -> float:
    p_s = compute_3to1_success_prob(f1, f2)
    if p_s <= 0:
        return f2
    f1_sq = f1 * f1
    one_minus_f1 = 1.0 - f1
    numerator = (
        f1_sq * f2
        + f1_sq * (1.0 - f2) / 27.0
        + 2.0 * f1 * one_minus_f1 * f2 / 9.0
        + one_minus_f1 * one_minus_f1 * f2 / 81.0
    )
    return min(1.0, max(0.0, numerator / p_s))


def compute_raw_channel_error(config: NoiseConfig) -> float:
    if config.channel_depolarization_error is not None and not isnan(config.channel_depolarization_error):
        return config.channel_depolarization_error
    p = config.physical_error
    return max(0.0, 1.0 - config.raw_epr_fidelity * (1.0 - p) * (1.0 - p))


def compute_remote_cnot_error_from_effective_channel(
    effective_channel_error: float,
    physical_error: float,
    accurate_rcx: bool,
) -> float:
    if not accurate_rcx:
        return max(0.0, effective_channel_error)
    return max(
        0.0,
        1.0 - (1.0 - effective_channel_error) * (1.0 - physical_error) * (1.0 - physical_error),
    )


def compute_distillation(protocol: str, rounds: int, config: NoiseConfig) -> DistillationResult:
    result = DistillationResult(
        protocol=protocol,
        output_fidelity=config.raw_epr_fidelity,
        success_probability=1.0,
        raw_pairs_consumed=1,
        distillation_time=0.0,
    )
    if protocol in ("none", "radar") or rounds == 0:
        return result

    f_raw = config.raw_epr_fidelity
    f_current = f_raw
    cumulative_success = 1.0
    total_pairs = 0

    for round_index in range(rounds):
        if protocol == "pumping_2to1":
            f_out = compute_2to1_fidelity(f_raw, f_current)
            p_s = compute_2to1_success_prob(f_raw, f_current)
            pairs_this_round = 2
        elif protocol == "pumping_3to1":
            f_out = compute_3to1_fidelity(f_raw, f_current)
            p_s = compute_3to1_success_prob(f_raw, f_current)
            pairs_this_round = 3
        elif protocol == "recurrence_2to1":
            f_out = compute_2to1_fidelity(f_current, f_current)
            p_s = compute_2to1_success_prob(f_current, f_current)
            pairs_this_round = 2 if round_index == 0 else (1 << (round_index + 1))
        elif protocol == "recurrence_3to1":
            f_out = compute_3to1_fidelity(f_current, f_current)
            p_s = compute_3to1_success_prob(f_current, f_current)
            pairs_this_round = 3 if round_index == 0 else (3 ** (round_index + 1))
        else:
            return result
        f_current = max(0.0, min(1.0, f_out))
        cumulative_success *= p_s
        total_pairs = pairs_this_round

    if protocol == "pumping_2to1":
        total_pairs = 2 * rounds
    elif protocol == "pumping_3to1":
        total_pairs = 3 * rounds

    result.output_fidelity = f_current
    result.success_probability = cumulative_success
    result.raw_pairs_consumed = max(1, total_pairs)
    if config.entanglement_rate > 0 and cumulative_success > 0:
        result.distillation_time = result.raw_pairs_consumed / (config.entanglement_rate * cumulative_success)
    return result


def compute_distillation_time(num_remote_cnots: int, distill: DistillationResult, config: NoiseConfig) -> float:
    if config.entanglement_rate <= 0:
        return 0.0
    expected_pairs = (num_remote_cnots * distill.raw_pairs_consumed) / distill.success_probability
    return expected_pairs / config.entanglement_rate


def check_timing_constraint(distillation_time: float, config: NoiseConfig) -> bool:
    return distillation_time <= config.measurement_delay


def evaluate_remote_cnot_noise_model(protocol: str, rounds: int, num_remote_cnots: int, config: NoiseConfig) -> SectionNoiseSummary:
    raw_channel_error = compute_raw_channel_error(config)
    distill = compute_distillation(protocol, rounds, config)
    distillation_time = compute_distillation_time(num_remote_cnots, distill, config)
    wants_distillation = protocol != "none" and rounds > 0
    branching = protocol_branching_factor(protocol)

    if wants_distillation:
        if is_recurrence_protocol(protocol):
            pairs_per_chain = branching ** rounds
        elif is_pumping_protocol(protocol):
            pairs_per_chain = branching
        else:
            pairs_per_chain = 1
        epr_slots = num_remote_cnots * pairs_per_chain * config.distillation_backup_batches
        distillation_qubits_required = 2 * epr_slots
        probability_all_fail = pow(1.0 - distill.success_probability, float(config.distillation_backup_batches))
        distilled_error = 1.0 - distill.output_fidelity
        effective_channel_error = (
            (1.0 - probability_all_fail) * distilled_error
            + probability_all_fail * raw_channel_error
        )
    else:
        epr_slots = 0
        distillation_qubits_required = 0
        probability_all_fail = 0.0
        distilled_error = 1.0 - distill.output_fidelity
        effective_channel_error = raw_channel_error

    remote_cnot_error = compute_remote_cnot_error_from_effective_channel(
        effective_channel_error,
        config.physical_error,
        config.accurate_rcx,
    )
    if config.local_baseline:
        remote_only_error = 0.0
        remote_cnot_error = 0.0
    elif config.monolithic_baseline:
        remote_only_error = 0.0
        remote_cnot_error = config.physical_error
    else:
        if config.interconnect_error > 0:
            remote_cnot_error = 1.0 - (1.0 - remote_cnot_error) * (1.0 - config.interconnect_error)
        remote_only_error = remote_cnot_error

    t_idle = max(config.measurement_delay, distillation_time)
    p_x = 0.0
    p_y = 0.0
    p_z = 0.0
    if config.t1_coherence_time > 0:
        p_x = (1.0 - exp(-t_idle / config.t1_coherence_time)) / 4.0
        p_y = p_x
    if config.t2_coherence_time > 0:
        p_z = (1.0 - exp(-t_idle / config.t2_coherence_time)) / 2.0 - p_x
        p_z = max(0.0, p_z)

    return SectionNoiseSummary(
        section_index=-1,
        repeat_count=0,
        remote_cnots_per_body=num_remote_cnots,
        required_epr_pairs_per_round=num_remote_cnots * max(1, distill.raw_pairs_consumed),
        effective_distillation_protocol=protocol,
        raw_channel_error=raw_channel_error,
        effective_channel_error=effective_channel_error,
        distilled_fidelity=distill.output_fidelity,
        distilled_error=distilled_error,
        distillation_success_probability=distill.success_probability,
        probability_all_distillation_fail=probability_all_fail,
        remote_cnot_error=remote_cnot_error,
        distillation_time_ns=distillation_time * 1e9,
        idling_time_us=t_idle * 1e6,
        p_x=p_x,
        p_y=p_y,
        p_z=p_z,
        timing_constraint_satisfied=check_timing_constraint(distillation_time, config),
        raw_pairs_per_distilled=distill.raw_pairs_consumed,
        distillation_backup_batches=config.distillation_backup_batches,
        distillation_epr_slots_required=epr_slots,
        distillation_qubits_required=distillation_qubits_required,
    )


def select_radar_distillation_protocol(rounds: int, num_remote_cnots: int, config: NoiseConfig) -> str:
    best_protocol = "none"
    best_score = float("inf")
    best_timing_safe = False
    best_pairs = 1 << 60
    best_time = float("inf")

    for candidate in ("none", "pumping_2to1", "pumping_3to1", "recurrence_2to1", "recurrence_3to1"):
        model = evaluate_remote_cnot_noise_model(candidate, rounds, num_remote_cnots, config)
        score = model.remote_cnot_error + model.p_x + model.p_y + model.p_z
        timing_safe = model.timing_constraint_satisfied
        better = score < best_score - 1e-15
        if not better and abs(score - best_score) <= 1e-15:
            if timing_safe != best_timing_safe:
                better = timing_safe
            elif model.raw_pairs_per_distilled != best_pairs:
                better = model.raw_pairs_per_distilled < best_pairs
            elif abs(model.distillation_time_ns - best_time) > 1e-9:
                better = model.distillation_time_ns < best_time
        if better:
            best_protocol = candidate
            best_score = score
            best_timing_safe = timing_safe
            best_pairs = model.raw_pairs_per_distilled
            best_time = model.distillation_time_ns
    return best_protocol


def effective_protocol(config: NoiseConfig, num_remote_cnots: int) -> str:
    if config.distillation_protocol != "radar":
        return config.distillation_protocol
    return select_radar_distillation_protocol(config.distillation_rounds, num_remote_cnots, config)


def infer_split_axis_and_value(circuit: stim.Circuit, requested_axis: str = "auto", requested_value: float | None = None) -> tuple[str, float]:
    coords = circuit.get_final_qubit_coordinates()
    if requested_axis not in ("auto", "x", "y"):
        raise ValueError(f"Unsupported split axis: {requested_axis}")

    if requested_axis == "x":
        axis = "x"
    elif requested_axis == "y":
        axis = "y"
    else:
        xs = sorted({coord[0] for coord in coords.values()})
        ys = sorted({coord[1] for coord in coords.values()})
        x_span = xs[-1] - xs[0]
        y_span = ys[-1] - ys[0]
        axis = "y" if y_span > x_span else "x"

    if requested_value is not None:
        return axis, requested_value

    values = sorted({coord[0] if axis == "x" else coord[1] for coord in coords.values()})
    return axis, values[len(values) // 2]


def _flat_remote_pairs(
    targets: Iterable[stim.GateTarget],
    coords: dict[int, list[float]],
    split_axis: str,
    split_value: float,
) -> list[int]:
    values = [target.value for target in targets]
    remote: list[int] = []
    axis_index = 0 if split_axis == "x" else 1
    for i in range(0, len(values), 2):
        a = values[i]
        b = values[i + 1]
        va = coords[a][axis_index]
        vb = coords[b][axis_index]
        crosses = (va <= split_value and vb > split_value) or (va > split_value and vb <= split_value)
        if crosses:
            remote.extend((a, b))
    return remote


def inject_network_noise(circuit: stim.Circuit, config: NoiseConfig) -> tuple[stim.Circuit, dict]:
    coords = circuit.get_final_qubit_coordinates()
    split_axis, split_value = infer_split_axis_and_value(
        circuit,
        requested_axis=config.split_axis,
        requested_value=config.split_x,
    )
    all_qubits = list(range(circuit.num_qubits))
    section_summaries: list[SectionNoiseSummary] = []
    has_repeat_sections = any(op.name == "REPEAT" for op in circuit)

    def transform_body(body: stim.Circuit, section_index: int, repeat_count: int) -> tuple[stim.Circuit, bool]:
        remote_cnots = 0
        for op in body:
            if op.name in REMOTE_2Q_GATES:
                remote_cnots += len(_flat_remote_pairs(op.targets_copy(), coords, split_axis, split_value)) // 2

        protocol = effective_protocol(config, remote_cnots)
        model = evaluate_remote_cnot_noise_model(protocol, config.distillation_rounds, remote_cnots, config)
        model.section_index = section_index
        model.repeat_count = repeat_count
        section_summaries.append(model)

        transformed = stim.Circuit()
        body_has_remote = remote_cnots > 0
        # Always inject idling at start of every repeat section (matches C++ which injects
        # PAULI_CHANNEL_1 into every REPEAT block, not just those with remote CNOTs).
        if model.p_x > 0 or model.p_y > 0 or model.p_z > 0:
            transformed.append("PAULI_CHANNEL_1", all_qubits, [model.p_x, model.p_y, model.p_z])

        for op in body:
            if op.name == "REPEAT":
                nested_body, _ = transform_body(op.body_copy(), len(section_summaries), op.repeat_count)
                transformed.append(stim.CircuitRepeatBlock(op.repeat_count, nested_body))
                continue
            transformed.append(op)
            if op.name in LOCAL_2Q_GATES and config.physical_error > 0 and not config.monolithic_baseline:
                # Local physical depolarizing noise on every 2-qubit gate (matches C++ inject_local_noise).
                all_targets = [t.value for t in op.targets_copy()]
                transformed.append("DEPOLARIZE2", all_targets, config.physical_error)
            if op.name in REMOTE_2Q_GATES:
                remote_pairs = _flat_remote_pairs(op.targets_copy(), coords, split_axis, split_value)
                if remote_pairs and model.remote_cnot_error > 0 and not config.monolithic_baseline:
                    transformed.append("DEPOLARIZE2", remote_pairs, model.remote_cnot_error)
        return transformed, body_has_remote

    def transform_flat_segment(segment: list[stim.CircuitInstruction], section_index: int) -> stim.Circuit:
        remote_cnots = 0
        for op in segment:
            if op.name in REMOTE_2Q_GATES:
                remote_cnots += len(_flat_remote_pairs(op.targets_copy(), coords, split_axis, split_value)) // 2

        protocol = effective_protocol(config, remote_cnots)
        model = evaluate_remote_cnot_noise_model(protocol, config.distillation_rounds, remote_cnots, config)
        model.section_index = section_index
        model.repeat_count = 1
        section_summaries.append(model)

        transformed = stim.Circuit()
        if model.p_x > 0 or model.p_y > 0 or model.p_z > 0:
            transformed.append("PAULI_CHANNEL_1", all_qubits, [model.p_x, model.p_y, model.p_z])

        for op in segment:
            transformed.append(op)
            if op.name in LOCAL_2Q_GATES and config.physical_error > 0 and not config.monolithic_baseline:
                all_targets = [t.value for t in op.targets_copy()]
                transformed.append("DEPOLARIZE2", all_targets, config.physical_error)
            if op.name in REMOTE_2Q_GATES:
                remote_pairs = _flat_remote_pairs(op.targets_copy(), coords, split_axis, split_value)
                if remote_pairs and model.remote_cnot_error > 0 and not config.monolithic_baseline:
                    transformed.append("DEPOLARIZE2", remote_pairs, model.remote_cnot_error)
        return transformed

    transformed = stim.Circuit()
    top_level_remote_pairs = 0
    repeat_section_index = 0
    if not has_repeat_sections:
        flat_segment: list[stim.CircuitInstruction] = []
        for op in circuit:
            if op.name == "RX" and flat_segment:
                transformed += transform_flat_segment(flat_segment, len(section_summaries))
                flat_segment = []
            flat_segment.append(op)
        if flat_segment:
            transformed += transform_flat_segment(flat_segment, len(section_summaries))
        summary = {
            "split_axis": split_axis,
            "split_value": split_value,
            "qubits": circuit.num_qubits,
            "detectors": circuit.num_detectors,
            "observables": circuit.num_observables,
            "top_level_remote_pairs": sum(section.remote_cnots_per_body for section in section_summaries),
            "config": asdict(config),
            "sections": [asdict(section) for section in section_summaries],
            "required_epr_pairs_per_round": max(
                (section.required_epr_pairs_per_round for section in section_summaries),
                default=0,
            ),
        }
        return transformed, summary

    for op in circuit:
        if op.name == "REPEAT":
            new_body, _ = transform_body(op.body_copy(), repeat_section_index, op.repeat_count)
            transformed.append(stim.CircuitRepeatBlock(op.repeat_count, new_body))
            repeat_section_index += 1
            continue
        transformed.append(op)
        if op.name in LOCAL_2Q_GATES and config.physical_error > 0 and not config.monolithic_baseline:
            all_targets = [t.value for t in op.targets_copy()]
            transformed.append("DEPOLARIZE2", all_targets, config.physical_error)
        if op.name in REMOTE_2Q_GATES:
            remote_pairs = _flat_remote_pairs(op.targets_copy(), coords, split_axis, split_value)
            top_level_remote_pairs += len(remote_pairs) // 2
            if remote_pairs:
                protocol = effective_protocol(config, len(remote_pairs) // 2)
                model = evaluate_remote_cnot_noise_model(protocol, config.distillation_rounds, len(remote_pairs) // 2, config)
                if model.remote_cnot_error > 0 and not config.monolithic_baseline:
                    transformed.append("DEPOLARIZE2", remote_pairs, model.remote_cnot_error)

    summary = {
        "split_axis": split_axis,
        "split_value": split_value,
        "qubits": circuit.num_qubits,
        "detectors": circuit.num_detectors,
        "observables": circuit.num_observables,
        "top_level_remote_pairs": top_level_remote_pairs,
        "config": asdict(config),
        "sections": [asdict(section) for section in section_summaries],
        "required_epr_pairs_per_round": (
            max((section.required_epr_pairs_per_round for section in section_summaries), default=top_level_remote_pairs)
        ),
    }
    return transformed, summary
