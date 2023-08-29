#!/usr/bin/env python
# Created on: Aug 27, 2023 at 1:28:20 PM
__author__ = "Michael Cuoco"

# configure logging
import logging

logger = logging.getLogger(__name__)

# Create a handler to output to the console
handler = logging.StreamHandler()

# Define the format
log_format = "[%(name)s: %(asctime)s] {%(thread)d} %(levelname)s - %(message)s"

# Configure the formatter with the provided format
formatter = logging.Formatter(log_format, datefmt="%m-%d %H:%M:%S")

# Attach the formatter to the handler
handler.setFormatter(formatter)

# Attach the handler to the logger
logger.addHandler(handler)

# Set the logging level
logger.setLevel(logging.INFO)
