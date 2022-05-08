"""
This script runs the ProPreFlask application using a development server.
"""

from os import environ
from ProPreFlask import app

if __name__ == '__main__':
    HOST = environ.get('SERVER_HOST', 'localhost')
    try:
        PORT = int(environ.get('SERVER_PORT', '56000'))
    except ValueError:
        PORT = 56000
    app.run(HOST, 57000)
