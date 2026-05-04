# BioNexus

## Setting up database

1) Copy `.env.example` -> `.env` and adjust

2) Start Postgres + Adminer:

    ```bash
    docker-compose --env-file .env.db -f docker-compose.db.yml up -d db adminer
    ```

    Adminer available at http://localhost:8080 (server: db, user: bionexus, db: bionexus)

3) Create/upgrade schema:

   ```bash
   set -a; source .env.db; set +a; python -m alembic -c alembic.ini upgrade head
   ```
   
   If you need to downgrade:
   ```bash
   set -a; source .env.db; set +a; python -m alembic -c alembic.ini downgrade -1
   ```

4) Load data with scripts in the `etl_scripts` folder.

5) Dump for deployment:

    Make sure the database is running in Docker, then:

    ```bash
    mkdir -p ~/Downloads/dumps && docker exec -i retromol-db-1 pg_dump -U bionexus -Fc -d bionexus -f /tmp/bionexus.dump && docker cp retromol-db-1:/tmp/bionexus.dump ~/Downloads/dumps/bionexus_$(date +%Y%m%d).dump
    ```