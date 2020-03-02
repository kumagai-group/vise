#!/usr/bin/env bash

if [ -z "$CI" ]; then
    echo "Will only continue on CI"
    exit
fi

if [[ $CIRCLE_BRANCH != "master" ]]; then
    echo "Will only continue for master builds"
    exit
fi

# build package and upload to private pypi index
touch ~/.pypirc
{
  echo "[distutils]"
  echo "index-servers = pypi-private"
  echo "[pypi-private]"
  echo "repository=https://$PYPI_HOST"
  echo "username=$PYPI_USERNAME"
  echo "password=$PYPI_PASSWORD"
} >> ~/.pypirc

python setup.py sdist upload -r pypi-private