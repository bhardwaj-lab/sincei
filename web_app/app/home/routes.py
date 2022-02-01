# -*- encoding: utf-8 -*-
"""
Copyright (c) 2019 - present AppSeed.us
"""

from app.home import blueprint
from flask import render_template, redirect, url_for, request, flash
from flask_login import login_required, current_user
from app import login_manager, db, Sample
from jinja2 import TemplateNotFound
import pickle

from bokeh.resources import INLINE

import sys
from .sincei_forms import form_scFilterStats, form_scPlotUMAP
from .sincei_functions import execute_command, fetch_results_scFilterStats, fetch_results_UMAP

# Helper - Extract current page name from request
def get_segment( request ):
    try:
        segment = request.path.split('/')[-1]
        if segment == '':
            segment = 'index'
        return segment
    except:
        return None

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
    form = form_scFilterStats(meta={'csrf': False})#, formdata=None
    flash(form.errors)
    if request.method == 'POST' and form.validate_on_submit():
        try:
            result = execute_command("scFilterStats", args=form)
            print(result)
            flash('Input submitted!!')
            #result=True
            if result:
                script, div = fetch_results_scFilterStats()
                return render_template(
                        'filter_stats.html',
                        form=form,
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
                form=form,
                plot_script="",
                plot_div="",
                js_resources=INLINE.render_js(),
                css_resources=INLINE.render_css(),
                ).encode(encoding='UTF-8')


@blueprint.route('/explore-output', methods = ['GET', 'POST'])
@login_required
def explore_umap():
    form = form_scPlotUMAP(meta={'csrf': False})#, formdata=None
    flash(form.errors)
    ## load the anndata object from output folder
    ad, df = load_anndata('filtered')
    if request.method == 'POST' and form.validate_on_submit():
        try:
            script, div = fetch_results_UMAP(ad, df)
            script_res, div_res = fetch_results_UMAP(ad, df, str(form.geneName.data))
            return render_template(
                        'explore-output.html',
                        form=form,
                        plot_script=script,
                        plot_div=div,
                        plot_script_res=script_res,
                        plot_div_res=div_res,
                        js_resources=INLINE.render_js(),
                        css_resources=INLINE.render_css(),
                        ).encode(encoding='UTF-8')
        except:
            return render_template('page-500.html'), 500
    else:
        script, div = fetch_results_UMAP()
        return render_template(
                    'explore-output.html',
                    form=form,
                    plot_script=script,
                    plot_div=div,
                    plot_script_res="",
                    plot_div_res="",
                    js_resources=INLINE.render_js(),
                    css_resources=INLINE.render_css(),
                    ).encode(encoding='UTF-8')
