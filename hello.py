from cloudant import Cloudant
from flask import Flask, render_template, request, jsonify
import atexit
import os
import json
import GDAL_Stuff

app = Flask(__name__, static_url_path='')

db_name = 'mydb'
client = None
db = None

if 'VCAP_SERVICES' in os.environ:
    vcap = json.loads(os.getenv('VCAP_SERVICES'))
    print('Found VCAP_SERVICES')
    if 'cloudantNoSQLDB' in vcap:
        creds = vcap['cloudantNoSQLDB'][0]['credentials']
        user = creds['username']
        password = creds['password']
        url = 'https://' + creds['host']
        client = Cloudant(user, password, url=url, connect=True)
        db = client.create_database(db_name, throw_on_exists=False)
elif "CLOUDANT_URL" in os.environ:
    client = Cloudant(os.environ['CLOUDANT_USERNAME'], os.environ['CLOUDANT_PASSWORD'], url=os.environ['CLOUDANT_URL'], connect=True)
    db = client.create_database(db_name, throw_on_exists=False)
elif os.path.isfile('vcap-local.json'):
    with open('vcap-local.json') as f:
        vcap = json.load(f)
        print('Found local VCAP_SERVICES')
        creds = vcap['services']['cloudantNoSQLDB'][0]['credentials']
        user = creds['username']
        password = creds['password']
        url = 'https://' + creds['host']
        client = Cloudant(user, password, url=url, connect=True)
        db = client.create_database(db_name, throw_on_exists=False)

# On IBM Cloud Cloud Foundry, get the port number from the environment variable PORT
# When running this app on the local machine, default the port to 8000
port = int(os.getenv('PORT', 8002))

@app.route('/Sweather')
def sweather():
  data = request.values
  if "long" in data.keys():
    if "lat" in data.keys():
      antiR = 10
      if "anticipWaterLevel" in data.keys():
        antiR = float(data["anticipWaterLevel"])
      print(data["long"],data["lat"])
      minRad = GDAL_Stuff.MIN_RADIUS
      if "minRadius" in data.keys():
        minRad = float(data["minRadius"])
      maxRad = GDAL_Stuff.MAX_RADIUS
      if "maxRadius" in data.keys():
        maxRad = float(data["maxRadius"])
      nReq = 1
      if "numToFind" in data.keys():
        nReq = int(data["numToFind"])
      divisions = GDAL_Stuff.TIFF_DIV
      if "threshDiv" in data.keys():
        threshDiv = int(data["threshDiv"])

      listGood = GDAL_Stuff.sweather_main(float(data["long"]),float(data["lat"]),anticipWater=antiR,minR=minRad,maxR=maxRad,nToSolve=nReq,threshDiv=threshDiv)
      print(listGood)
  return jsonify(listGood)
  

@app.route('/')
def root():    
    data = request.data
    print(data)
    return app.send_static_file("static/index.html")

# /* Endpoint to greet and add a new visitor to database.
# * Send a POST request to localhost:8000/api/visitors with body
# * {
# *     "name": "Bob"
# * }
# */
@app.route('/api/visitors', methods=['GET'])
def get_visitor():
    if client:
        return jsonify(list(map(lambda doc: doc['name'], db)))
    else:
        print('No database')
        return jsonify([])

# /**
#  * Endpoint to get a JSON array of all the visitors in the database
#  * REST API example:
#  * <code>
#  * GET http://localhost:8000/api/visitors
#  * </code>
#  *
#  * Response:
#  * [ "Bob", "Jane" ]
#  * @return An array of all the visitor names
#  */
@app.route('/api/visitors', methods=['POST'])
def put_visitor():
    user = request.json['name']
    data = {'name':user}
    if client:
        my_document = db.create_document(data)
        data['_id'] = my_document['_id']
        return jsonify(data)
    else:
        print('No database')
        return jsonify(data)

@atexit.register
def shutdown():
    if client:
        client.disconnect()

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=port, debug=True)
    