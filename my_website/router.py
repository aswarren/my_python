from flask import Flask
from flask import render_template
from flask import request
app = Flask(__name__)
import re
import protein_to_gene_ids

@app.route('/')
def hello_world():
    return 'Hello World!'

@app.route('/tools/geneid/', methods=['POST', 'GET'])
def get_geneids():
    if request.method == 'POST':
         if re.match("[^@]+@[^@]+\.[^@]+", request.form['email']) and len(request.form['id_text']):
             return geneid_results(request.form['email'], request.form['id_text'])
    return render_template('geneid.html')

def geneid_results(email, id_text):
    results=protein_to_gene_ids.process_text(email, id_text)
    return render_template('geneid.html', results=results)

if __name__ == '__main__':
    #app.debug=True
    app.run(host='0.0.0.0')
