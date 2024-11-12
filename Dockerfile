FROM python:3.9-slim

WORKDIR /app

COPY requirements.txt .

# Install system dependencies required for RDKit and other chemical tools
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libxext6 \
    libgl1-mesa-glx \
    && rm -rf /var/lib/apt/lists/*

# If you're using OpenBabel, you might also need:
RUN apt-get update && apt-get install -y \
    openbabel \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "80"]