# Ultra-light image for maintenance loop
FROM python:3.10-slim

WORKDIR /app

# Unbuffered logs so they appear in docker logs immediately
ENV PYTHONUNBUFFERED=1

# Copy only server code (includes maintenance.py and routes/)
COPY src/server /app

# Install the same dependencies as the backend (could be slimmed down)
RUN pip install --no-cache-dir -r requirements.maintenance.txt

# Run the maintenance loop
CMD ["python", "-m", "maintenance"]