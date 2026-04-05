-- Create a dedicated read-only role for the app
DO $$
BEGIN
  IF NOT EXISTS (SELECT FROM pg_roles WHERE rolname = 'app_ro') THEN
    CREATE ROLE app_ro LOGIN PASSWORD 'apppass_ro';
  ELSE
    ALTER ROLE app_ro WITH PASSWORD 'apppass_ro';
  END IF;
END$$;

-- Harden the role so transactions default to read-only
ALTER ROLE app_ro SET default_transaction_read_only = on;

-- Ensure ownership stays with the superuser created by POSTGRES_USER
-- Revoke public writes; grant minimal read permissions
REVOKE CREATE ON SCHEMA public FROM PUBLIC;
GRANT USAGE ON SCHEMA public TO app_ro;

-- Grant connect on the current database
DO $$
BEGIN
  EXECUTE format('GRANT CONNECT ON DATABASE %I TO app_ro;', current_database());
END $$;

-- Grant SELECT on all existing tables
GRANT SELECT ON ALL TABLES IN SCHEMA public TO app_ro;
GRANT SELECT ON ALL SEQUENCES IN SCHEMA public TO app_ro;

-- Ensure future objects are readable but not writable
ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT SELECT ON TABLES TO app_ro;
ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT SELECT ON SEQUENCES TO app_ro;