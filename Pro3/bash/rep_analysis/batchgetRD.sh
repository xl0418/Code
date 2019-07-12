#!/bin/bash

bash batchuni2win.sh &&
bash batchextractallDs.sh &&
bash batchextractallRs.sh &&
bash batchDlast_csv.bash &&
bash batchRlast_csv.bash &&
bash batchloopextractEve_Turn_csv.bash &&
bash batchcleanfiles.bash
