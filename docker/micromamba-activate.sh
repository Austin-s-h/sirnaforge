# Nextflow sources this file with the requested process environment as its
# argument. Define micromamba's shell integration in the current bash process,
# then activate that environment without invoking a subprocess shell.
eval "$(/usr/local/bin/micromamba shell hook --shell bash)"
micromamba activate "$@"
