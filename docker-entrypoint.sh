#!/bin/bash
set -e

rm -f ~/.julia/compiled/v1.12/PeriLab/*.ji

exec "$@"
