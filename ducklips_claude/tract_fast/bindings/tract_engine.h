/* The in-process engine: the sustain tool's run (tract::SustainRun) behind a
 * plain C interface for ctypes. A render's inputs are given as numbers and
 * arrays (never text to parse); its per-step outputs are written into the
 * caller's buffers. A run can be copied at any step (a snapshot) and the
 * copy resumed, with inputs that act from that step on changed.
 *
 * Settings are named by the sustain tool's environment variables (VISCOUS,
 * HOPLOSS, SEAL, REVJET, ...), so the mapping is the tool's own; arrays by
 * the files those variables name (AREAFILE, AREASEQ, VELUMSEQ, FOLDSEQ,
 * FOLDKEYS, JETSEQ, FORCEWAVE, FORCEFLOW, AREAFILE2).
 *
 * Every function returning int returns 0 on success, nonzero on an unknown
 * name or bad shape (tract_error then describes it). */
#ifndef TRACT_ENGINE_H
#define TRACT_ENGINE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct tract_settings tract_settings;
typedef struct tract_run tract_run;

/* The last error's text (thread-local), empty when none. */
const char* tract_error(void);

tract_settings* tract_settings_new(void);
tract_settings* tract_settings_copy(const tract_settings* settings);
void tract_settings_free(tract_settings* settings);

/* The command line's arguments: PITCH, LUNG, SECONDS. And any numeric
 * environment variable (PITCHEND RELEASE ASPIRATION FRICATION FLUX CLOSING
 * TREMOR ADDTREMOR BREATHTREMOR PITCHFB PORTRESIST REVJET REVJETCLOSE
 * FOLDTURB SKIN WALLS WALLMASS AXIAL WALLVAR VISCOUS ACTIVELENGTH STATICAREA
 * REVNOSE SEAL NOSEVISCOUS NOSESHELF NOSEHOP HOPLOSS CHINK VELUM OBSERVE
 * RAWDOUBLE RAWCORD): the value the tool would read from its text. */
int tract_settings_number(tract_settings* settings, const char* name, double value);
/* "first,second" variables: NOSERAD REST FOLDDAMP GLIDE. */
int tract_settings_pair(tract_settings* settings, const char* name, double first, double second);
/* Text variables: BRANCHES, and the output paths WAVE (the tool's OUT.wav;
 * unset: none written) SOURCE NOSEDUMP AREAS WALLDUMP RAW PORTWAVES PORTDUMP
 * GEODUMP (written exactly as the tool writes them). */
int tract_settings_text(tract_settings* settings, const char* name, const char* text);
/* Arrays (row-major, rows x columns), as the tool reads its files:
 *   AREAFILE  1 x (1 + n): the section length, then n areas
 *   AREAFILE2 1 x (1 + n): (the length is ignored, as the tool ignores it)
 *   AREASEQ   rows x (1 + n): time, then n areas
 *   VELUMSEQ, JETSEQ  rows x 2;  FOLDKEYS rows x 4;  FOLDSEQ rows x 5
 *   FORCEWAVE, FORCEFLOW  1 x steps (NaN where none) */
int tract_settings_array(tract_settings* settings, const char* name, const double* values, int64_t rows,
                         int64_t columns);
/* Variables as text, parsed exactly as the tool parses its environment
 * (files they name are read), for the variables a caller's environment
 * passes through. Apply before the numbers and arrays the caller sets
 * itself: each set variable replaces what is there. AREAFILE first. */
int tract_settings_environment(tract_settings* settings, int count, const char* const* names,
                               const char* const* values);
/* Steps a render of these settings takes (seconds * 96000, truncated). */
int64_t tract_settings_steps(const tract_settings* settings);
/* The vocal tube's node count for these settings (GEODUMP's rows per step). */
int64_t tract_settings_vocal_nodes(const tract_settings* settings);

/* The first step whose inputs differ between two settings (the controls
 * applied before it steps, its forced wave and observed flow): a run of
 * `second` may resume from any state a run of `first` passed through at or
 * before that step, and its outputs before it are `first`'s. -1 when the
 * two build different bodies (or a dump written after the first step is
 * set): not even the state before the first step is shared. */
int64_t tract_settings_divergence(const tract_settings* first, const tract_settings* second);

tract_run* tract_run_new(const tract_settings* settings);
tract_run* tract_run_copy(const tract_run* run);
void tract_run_free(tract_run* run);
int64_t tract_run_step(const tract_run* run);
int64_t tract_run_total(const tract_run* run);  /* as of its last advance */
/* Steps the run to `end` (at most the steps of `settings`) with these
 * settings' inputs (the run's own or changed ones that act from its step
 * on; the run's length and RELEASE become theirs). Per step from
 * the run's step, rows are written from the start of each buffer given
 * (NULL: not wanted): raw 2 doubles (the glottal flow entering and the
 * radiated flow, RAWDOUBLE's pair), port_waves 6 (PORTWAVES' row),
 * geometry vocal_nodes x 12 (GEODUMP's rows). Returns the steps written
 * (with a non-finite output the step that made it is the last, and
 * *status is 1; else 0). */
int64_t tract_run_advance(tract_run* run, const tract_settings* settings, int64_t end, double* raw,
                          double* port_waves, double* geometry, int32_t* status);
/* The tool's ending: the WAVE (and SOURCE) written when their paths are set. */
void tract_run_finish(const tract_run* run, const tract_settings* settings);
/* The nose as NOSEDUMP gives it, for the run as it stands: header 4 values
 * (port node, node count, hop low-pass, viscous scale), then per nose node
 * 12 (area, 0, link resistance, the lining's nine). Returns the nose's node
 * count; writes at most `capacity` nodes. */
int64_t tract_run_nose(const tract_run* run, double* header, double* rows, int64_t capacity);

/* One whole render (new, advance to the end, finish, free). */
int64_t tract_render(const tract_settings* settings, double* raw, double* port_waves, double* geometry,
                     int32_t* status);

#ifdef __cplusplus
}
#endif

#endif
