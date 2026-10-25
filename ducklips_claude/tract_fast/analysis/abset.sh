#!/bin/bash
# Sentence 1 under each suspect removed, for the owner's ear.
cd /tmp/es; S=~/ducklips/tract/build/sheet_speak; mkdir -p ab; rm -f ab/*.wav
r() { $S v1u.score ab/$1.wav ${@:2} >/dev/null 2>&1; }
r A_current
r B_walls_silent --skin 0
r C_stiff_walls --tract-mass 10 --skin 0
r D_no_muscle_wander --hold-noise 0 --hold-common 0 --tissue 0
r E_softer_breath --lung-scale 0.6
r F_no_pitch_selfcorrection --pitch-feedback 0
r G_no_lung_branch --lung-resistance 0
r H_old_airway --subglottal-narrowing 1 --lung-scale 1
$S v1flat.score ab/I_flat_pitch_command.wav >/dev/null 2>&1
r J_everything_off --skin 0 --tract-mass 10 --hold-noise 0 --hold-common 0 --tissue 0 --pitch-feedback 0 --elongation 0 --devoicing-raise 0
ls ab
