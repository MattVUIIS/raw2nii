#!/bin/bash
zip -r raw2nii __main__.py project
echo '#!/usr/bin/env python' | cat - raw2nii.zip > raw2nii
chmod u+x raw2nii
rm raw2nii.zip
