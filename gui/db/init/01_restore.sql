#!/bin/sh
# Restores a Postgres database from a .dump file mounted at /docker-entrypoint-initdb.d/dump.dump
# Make sure it is executable: chmod +x db/init/01_restore.sh
set -eu

FILE="/docker-entrypoint-initdb.d/dump.dump"

if [ ! -f "$FILE" ]; then
  echo "No dump file found at $FILE"
  exit 1
fi

case "$FILE" in
  *.dump)
    echo "Restoring pg_dump archive: $FILE"
    pg_restore -U "$POSTGRES_USER" -d "$POSTGRES_DB" -v "$FILE"
    echo "Restore complete."
    ;;
  *)
    echo "Unsupported dump file type (expected .dump): $FILE"
    exit 1
    ;;
esac