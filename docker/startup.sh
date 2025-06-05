#!/bin/bash

# Pull or clone bipm_g2o from the master branch
if [ ! -d /workspace/ecg2o ]; then
    git clone -b master https://github.com/snt-arg/ecg2o.git
else
    echo "Repo already exists, updating..."
    cd /workspace/ecg2o
    git checkout master
    git pull origin master
fi

# Enter the project directory
cd /workspace/ecg2o

# Start an interactive shell
exec /bin/bash

