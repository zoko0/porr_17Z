#!/bin/bash
chpl -o execute test.chpl
./execute
rm execute
