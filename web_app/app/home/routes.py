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

# for bokeh plotting
from bokeh.models import ColumnDataSource, Div, Select, Slider, TextInput
from bokeh.io import curdoc
from bokeh.resources import INLINE
from bokeh.embed import components
from bokeh.plotting import figure, output_file, show
import pandas as pd


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
        print("not found! "+template)
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

@blueprint.route('/filter_stats', methods = ['GET', 'POST'])
@login_required
def show_plots():
    if request.method == 'POST':
        try:
            def getFilterStats(txtpath="/Users/vivek/programs/sincei/web_app/example_data/scFilterStats.txt"):
                res = pd.read_csv(txtpath, sep="\t", index_col=0)
                return res

            df = getFilterStats()
            df['color']="#000000"
            pretty_labels={}
            for x in df.columns:
                pretty_labels[" ".join(x.split("_"))] = x

            source = ColumnDataSource(df)
            xlabel="Total sampled"
            ylabel="Wrong motif"
            fig = figure(plot_height=500, plot_width=500, tooltips=[("ID", "@index")])# level_0 refers to "index"
            fig.circle(x=pretty_labels[xlabel], y=pretty_labels[ylabel], source=source, size=8, color="color", line_color=None)
            fig.xaxis.axis_label = xlabel
            fig.yaxis.axis_label = ylabel
            script, div = components(fig)
            return render_template(
                    'filter_stats.html',
                    plot_script=script,
                    plot_div=div,
                    js_resources=INLINE.render_js(),
                    css_resources=INLINE.render_css(),
                    ).encode(encoding='UTF-8')
        except:
            return render_template('page-500.html'), 500
    else:
        return render_template(
                'filter_stats.html',
                plot_script="",
                plot_div="",
                js_resources=INLINE.render_js(),
                css_resources=INLINE.render_css(),
                ).encode(encoding='UTF-8')

# Helper - Extract current page name from request
def get_segment( request ):
    try:
        segment = request.path.split('/')[-1]
        if segment == '':
            segment = 'index'
        return segment
    except:
        return None
