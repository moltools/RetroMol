# Python backend with gunicorn
FROM condaforge/mambaforge:24.9.2-0

# Let us use bash -lc in RUN commands (so conda works nicely)
SHELL ["/bin/bash", "-lc"]

# Create non-root user for saner permissions
ARG USERNAME=app
ARG USER_UID=1000
ARG USER_GID=$USER_UID
RUN groupadd --gid $USER_GID $USERNAME \
 && useradd  --uid $USER_UID --gid $USER_GID -m $USERNAME

 # Switch to /app as working dir
WORKDIR /app

# System deps (git; libxrender1/libxext6/libsm6 are runtime deps of rdkit.Chem.Draw,
# which dlopens X11 libs for 2D rendering even though the container is headless)
RUN apt-get update && apt-get install -y --no-install-recommends build-essential git libxrender1 libxext6 libsm6 && rm -rf /var/lib/apt/lists/*

# Copy env + requirements before env creation for caching
COPY gui/src/server/environment.backend.yml /app/
COPY gui/src/server/requirements.backend.txt /app/

# Build conda env (includes pip deps via environment.backend.yml)
RUN mamba env create -n retromol-gui -f environment.backend.yml && conda clean -afy

# copy retromol package source + pyproject
COPY pyproject.toml /app/pyproject.toml
COPY README.md /app/README.md
COPY src /app/src

# install local retromol into the conda env
RUN conda run -n retromol-gui pip install -e /app

# Copy backend; Flask code lives in src/server
# Set ownership to non-root user
COPY gui/src/server /app

# Ensure model/cache dirs exist
RUN mkdir -p /app/models /app/cache && chown -R ${USER_UID}:${USER_GID} /app

# Runtime env vars
ENV CACHE_DIR=/app/cache \
    PYTHONUNBUFFERED=1 \
    LOG_LEVEL=INFO \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 \
    JOBLIB_MMAP_MODE=r \
    PORT=4000

# gunicorn entry (Flask app is created in routes/app.py as `app`)
# expose 4000 internally on the network
EXPOSE 4000

# User the unprivileged user at runtime
USER $USERNAME

# Run everything inside the conda env without manual activation
ENTRYPOINT ["conda", "run", "-n", "retromol-gui", "--no-capture-output"]

# Worker/thread/timeout sizing lives in gunicorn.conf.py (env-var overridable via
# gui/docker/backend.env), not hardcoded here.
# Let Flask/gunicorn find the app: "app:app"
CMD ["gunicorn", "-c", "gunicorn.conf.py", "app:app"]

# /api/ready checks both DuckDB and Redis connectivity -- a more meaningful signal
# than "the process is up" for whether this container should receive traffic. Uses
# stdlib urllib rather than curl/wget (neither is installed in this image, and
# adding one just for a healthcheck isn't worth the extra apt layer).
HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
    CMD conda run -n retromol-gui --no-capture-output python -c \
    "import urllib.request; urllib.request.urlopen('http://localhost:4000/api/ready', timeout=3)" || exit 1