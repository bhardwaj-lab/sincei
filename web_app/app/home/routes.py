# -*- encoding: utf-8 -*-
"""
Copyright (c) 2019 - present AppSeed.us
"""

from app.home import blueprint
from flask import render_template, redirect, url_for, request
from flask_login import login_required, current_user
from app import login_manager, db, Sample
from jinja2 import TemplateNotFound
import pickle

@blueprint.route('/index')
@login_required
def index():

    return render_template('index.html', segment='index')

@blueprint.route('/<template>')
@login_required
def route_template(template):
    try:
        if not template.endswith( '.html' ):
            template += '.html'
        # Detect the current page
        segment = get_segment( request )
        # Serve the file (if exists) from app/templates/FILE.html
        return render_template( template, segment=segment )
    except TemplateNotFound:
        return render_template('page-404.html'), 404
    except:
        return render_template('page-500.html'), 500

@blueprint.route('/uploader', methods = ['GET', 'POST'])
@login_required
def uploader_file():
   if request.method == 'POST':
      try:
          f = request.files.getlist("file[]")
          for name in f:
              filename = name.filename
              samplename = name.filename.split('.')[0]
              # insert data into the 'sample' SQLite table
              samples = Sample(filename, samplename)
              # Flask-SQLAlchemy magic adds record to database
              db.session.add(samples)
          db.session.commit()
          samples_out = Sample.query.order_by(Sample.samplename).all()
          return render_template('file_upload.html', samples=samples_out)
      except:
          return render_template('page-500.html'), 500

# Helper - Extract current page name from request
def get_segment( request ):
    try:
        segment = request.path.split('/')[-1]
        if segment == '':
            segment = 'index'
        return segment
    except:
        return None
