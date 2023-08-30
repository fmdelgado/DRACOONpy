# Use an official Python runtime as a parent image
FROM python:3.10-slim-buster

RUN apt update && apt dist-upgrade -y && apt install -y build-essential

# Set the working directory in the container to /app
WORKDIR /app

# Add the current directory contents into the container at /app
ADD app/ .

RUN pip install --upgrade pip

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Make port 8031 available to the world outside this container
EXPOSE 8031

# Run app.py when the container launches
CMD streamlit run dracoon_app.py --browser.serverAddress prototypes.cosy.bio --server.enableCORS=false --server.enableXsrfProtection=false --server.enableWebsocketCompression=false
