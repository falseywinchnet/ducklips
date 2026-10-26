/* The tube peeled from the sound, time-major (NOTES_residual_speed.md):
 * observer.StatesPeel's and reverse_core.SidePeel's problems, the same
 * values to the bit, run as a skewed wavefront over all layers at once
 * instead of node by node in chunks.
 *
 * A layer j (0 <= j <= n - 2) is link j and node j. Link j at step t reads
 * the wave arriving from node j + 1 at t and t - 1 and node j + 1's
 * pressure at t - 1; node j at t reads link j at t and t + 1. Layer j's
 * link therefore runs d_j = 1 + 2 (n - 2 - j) wavefronts behind the open
 * end (node n - 1, at wavefront W step W), and node j one step behind its
 * link: within one wavefront every layer reads only what earlier
 * wavefronts wrote, so the layers are independent (the processor overlaps
 * their division chains, which in the node-major peel run in series).
 *
 * Every value is computed by the numba kernels' expression on the same
 * inputs (this file is compiled with -ffp-contract=off: no multiply-add is
 * fused, as numba fuses none). */
#ifndef PEEL_CORE_H
#define PEEL_CORE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Which numba kernel's expressions: 0 observer._states_link (the shelf,
 * the record, a start from a forward state), 1 reverse_core._side_link
 * (no side port, no stop node). */
enum { PEEL_VARIANT_STATES = 0, PEEL_VARIANT_SIDE = 1 };

/* The peel's fixed inputs. Borrowed for the duration of each call. */
typedef struct peel_inputs {
  const double* geometry;     /* steps x nodes x 12 (moving_geometry's GEODUMP rows), C order */
  const double* sound;        /* steps: the flow radiated */
  int64_t steps;
  int64_t nodes;              /* n >= 3 */
  int64_t t_start;            /* first step peeled (STATES; 0 for SIDE) */
  int64_t record;             /* the forward state after this step into `recorded` (-1: none) */
  const double* init;         /* nodes x 19: the forward state after t_start - 1 (STATES with t_start > 0) */
  double* recorded;           /* nodes x 19 (STATES with record >= 0) */
  double* flow_out;           /* steps: node 0's flow entering (written where node 0 runs) */
  double* pressure_out;       /* steps: node 0's pressure */
  double hop;
  double hop_zero;
  double end_extra;
  double jet_gain;            /* 0: no reversible jet */
  double jet_closing;
  double jet_bound;
  const double* gate;         /* gate_length values (has_gate) */
  const double* noise;        /* noise_links x noise_steps (kernels_jet.noise_table), C order */
  int64_t gate_length;
  int64_t noise_links;
  int64_t noise_steps;
  int64_t jet_seed;
  int32_t kinetic;
  int32_t has_gate;
  int32_t variant;            /* PEEL_VARIANT_* */
  int32_t reserved;
} peel_inputs;

/* Doubles of a peel's state (layers, linings, open end, wave rings). */
int64_t peel_state_size(int64_t nodes);

/* The state before the first wavefront (zeros, or `init`'s forward state
 * at t_start for STATES; the jet's coefficient caches empty). */
void peel_begin(const peel_inputs* inputs, double* state);

/* Wavefronts begin..target - 1 (wavefront W runs the open end at step W).
 * With snap_every > 0, the state at each wavefront W in (begin, target]
 * that is a multiple of snap_every (and below peel_complete) is copied to
 * snapshots + (W / snap_every) * peel_state_size (capacity snap_slots
 * slots). Inputs read: steps below target only. */
void peel_advance(const peel_inputs* inputs, double* state, int64_t begin, int64_t target, double* snapshots,
                  int64_t snap_every, int64_t snap_slots);

/* The wavefront at which every node's every step is done. */
int64_t peel_complete(int64_t steps, int64_t nodes);

/* The moving geometry: reverse_core._geometry_steps to the bit (the GEODUMP
 * rows of steps begin..end - 1 into `out`, the recurrence carried in `state`
 * and `counters`, its state before every step that is a multiple of
 * snap_every into snaps / snap_counters), and per step what the peel's
 * callers would read from the rows in further passes: the least area
 * (observer.least_area), the first node's area, and the reversible jet's
 * narrowest node and its area (reverse_core._narrowest of the step's own
 * areas). Borrowed for the call. */
typedef struct geometry_inputs {
  const double* rest;         /* ceil(steps / rest_every) x nodes: the rest areas set before each step */
  const double* materials;    /* nodes x 9 (reverse_core.material_table) */
  double* out;                /* steps x nodes x 12 */
  double* state;              /* nodes x 10 (reverse_core.GEOMETRY_STATE) */
  double* counters;           /* 2: the ramp's steps left, the glide's countdown */
  double* snaps;              /* slots x nodes x 10 (snap_every > 0) */
  double* snap_counters;      /* slots x 2 */
  double* least;              /* steps */
  double* first_area;         /* steps */
  int64_t* jet_node;          /* steps */
  double* jet_area;           /* steps */
  int64_t steps;
  int64_t nodes;
  int64_t begin;
  int64_t end;
  int64_t snap_every;
  int64_t rest_every;         /* steps a row of `rest` holds for (1: a row a step; the command's tick: 96) */
  double viscous;
} geometry_inputs;

void geometry_steps(const geometry_inputs* inputs);

/* The tube run forward over steps t0..t1 - 1 (observer._forward_steps, a
 * hidden window's run: the folds' flow in, the lips' flow and the first
 * node's pressure out), from and into `state` (nodes x 19, observer's
 * state rows), to the bit. Borrowed for the call. */
typedef struct forward_inputs {
  const double* geometry;     /* steps x nodes x 12 */
  const double* flow_in;      /* steps */
  double* state;              /* nodes x 19 */
  double* lips_out;           /* steps */
  double* port_pressure;      /* steps */
  const double* gate;         /* gate_length (has_gate) */
  const double* noise;        /* noise_links x noise_steps */
  int64_t steps;
  int64_t nodes;
  int64_t t0;
  int64_t t1;
  int64_t gate_length;
  int64_t noise_links;
  int64_t noise_steps;
  int64_t jet_seed;
  double hop;
  double hop_zero;
  double end_extra;
  double seal;
  double jet_gain;
  double jet_closing;
  double jet_bound;
  int32_t kinetic;
  int32_t seal_blend;
  int32_t has_gate;
  int32_t reserved;
} forward_inputs;

void forward_steps(const forward_inputs* inputs);

#ifdef __cplusplus
}
#endif

#endif
