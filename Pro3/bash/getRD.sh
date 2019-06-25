#!/bin/bash

bash uni2win.sh &&
bash extractallDs.sh &&
bash extractallRs.sh &&
bash Dlast_csv.bash &&
bash Rlast_csv.bash &&
bash loopextractEve_Turn_csv.bash &&
bash cleanfiles.bash
