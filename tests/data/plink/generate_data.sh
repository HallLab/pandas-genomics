#!/bin/bash

# simulate with plink
plink --simulate medium.sim --simulate-ncases 100 --simulate-ncontrols 500 --make-bed --out plink_test_medium
plink --simulate small.sim --simulate-ncases 50 --simulate-ncontrols 100 --make-bed --out plink_test_small