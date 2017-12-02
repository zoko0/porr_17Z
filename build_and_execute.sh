#!/bin/bash
sudo docker run --rm -v "$PWD":/usr/src/myapp -w /usr/src/myapp chapel/chapel chpl -o execute main_loop.chpl
sudo docker run --rm -v "$PWD":/usr/src/myapp -w /usr/src/myapp chapel/chapel ./execute
sudo rm execute
