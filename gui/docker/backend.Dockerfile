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

# System deps (psycopg binary + git)
RUN apt-get update && apt-get install -y --no-install-recommends build-essential libpq-dev git && rm -rf /var/lib/apt/lists/*

# Copy env + requirements before env creation for caching
COPY src/server/environment.backend.yml /app/
COPY src/server/requirements.backend.txt /app/

# Build conda env (includes pip deps via environment.backend.yml)
RUN mamba env create -n retromol-gui -f environment.backend.yml && conda clean -afy

# Copy backend; Flask code lives in src/server
# Set ownership to non-root user
COPY src/server /app

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

# Let Flask/gunicorn find the app: "app:app"
CMD ["gunicorn", "-w", "1", "--threads", "4", "-b", "0.0.0.0:4000", "--access-logfile", "-", "--error-logfile", "-", "--log-level", "info", "--timeout", "120", "app:app"]