#!/usr/bin/env python3

from __future__ import print_function
from base64 import b64encode
import json
from urllib import request, parse

def get_access_token(client_id, client_secret, idt_username, idt_password):
    """
    Create the HTTP request, transmit it, and then parse the response for the 
    access token.
    
    The body_dict will also contain the fields "expires_in" that provides the 
    time window the token is valid for (in seconds) and "token_type".
    """

    # Construct the HTTP request
    authorization_string = b64encode(bytes(client_id + ":" + client_secret, "utf-8")).decode()
    request_headers = { "Content-Type" : "application/x-www-form-urlencoded",
                        "Authorization" : "Basic " + authorization_string }
                    
    data_dict = {   "grant_type" : "password",
                    "scope" : "test",
                    "username" : idt_username,
                    "password" : idt_password }
    request_data = parse.urlencode(data_dict).encode()

    post_request = request.Request("https://eu.idtdna.com/Identityserver/connect/token", 
                                    data = request_data, 
                                    headers = request_headers,
                                    method = "POST")

    # Transmit the HTTP request and get HTTP response
    response = request.urlopen(post_request)

    # Process the HTTP response for the desired data
    body = response.read().decode()
    
    # Error and return the response from the endpoint if there was a problem
    if (response.status != 200):
        raise RuntimeError("Request failed with error code:" + response.status + "\nBody:\n" + body)
    
    body_dict = json.loads(body)
    return body_dict["access_token"]
    
    
if __name__ == "__main__":
    client_id = "Arne"
    client_secret = "38a9cb9a-58a2-45e5-b7f9-b2a5d878d66a"
    idt_username = "ArBlom"
    idt_password = "51kCH8KSvD1F5*;"
    
    token = get_access_token(client_id, client_secret, idt_username, idt_password)
    print(token)


import requests

# Your API token (replace 'YOUR_API_TOKEN' with your actual token)
api_token = token

# Base URL for the IDT API (this is an example; you should check the actual base URL from IDT's documentation)
base_url = 'https://eu.idtdna.com/restapi'

# Specific endpoint for the OligoAnalyzer (adjust according to the actual endpoint provided by IDT)
endpoint = '/v1/OligoAnalyzer/TmMisMatch'

# Full URL
url = f"{base_url}{endpoint}"

# Headers including the Authorization token
headers = {
    'Authorization': f'Bearer {api_token}',
    'Content-Type': 'application/json'
}

# Parameters for the OligoAnalyzer (example, you should customize based on your needs)
params = {
  "Settings": {  "Sequence": "GTCTTTACTCACCTGTAGATGCCT",
    "NaConc": 50,
    "MgConc": 3,
    "dNTPsConc": 0.25,
    "OligoConc": 1.2,
    "SequenceType": "DNA"
  },
  "Sequence2": "GTATTTACTCACCTGTAGATGCCT"
}

# Send a POST request (or GET, depending on the API specification)
response = requests.post(url, headers=headers, json=params)

# Check the response status
if response.status_code == 200:
    # Parse the JSON response
    data = response.json()
    print(data)
else:
    print(f"Error: {response.status_code}")
    print(response.text)