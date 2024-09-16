# Use an official Python runtime as a parent image
FROM python:3.11.4


# Set the working directory in the container
WORKDIR /Users/utkarshprajapati/Desktop/Columbia Course Files/RA - Climate Modelling/climate_code_conversion/_codegen_claude


# Copy the requirements file into the container at /usr/src/app
COPY requirements.txt ./


# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt


# Copy the rest of the application into the container
COPY . .

# To check the installed packages
RUN pip install scipy
