#!/bin/sh -eux
useradd --no-create-home --uid 1000 dockuser
mkdir -p /workspace
chown -R dockuser:dockuser /workspace /autodocker
rm -f /setup.sh
