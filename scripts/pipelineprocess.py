from FlaskApp import app, db
from flask import flash
import os
import subprocess
from time import sleep

def clean_gwas(filename,code,form):
    flash('Cleaning file...')
    statinfo = os.stat(os.path.join(app.config['UPLOAD_FOLDER'], filename))
    if statinfo.st_size >= 1048576:
        cmd = [os.path.join(app.config['PYTHON_FOLDER'],'python'),
               os.path.join(app.config['SCRIPT_FOLDER'],'cleanSumstats.py'), 
               os.path.join(app.config['UPLOAD_FOLDER'], filename),
               os.path.join(app.config['CLEAN_FOLDER'], code)]
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        (output, err) = process.communicate()
        p_status = process.wait()
        statinfo = os.stat(os.path.join(app.config['CLEAN_FOLDER'], code+'.gz')) 
#        if True:
        if statinfo.st_size >= 1048576:
            cmd = [os.path.join(app.config['PYTHON_FOLDER'],'python'),
                   os.path.join(app.config['SCRIPT_FOLDER'],'completeclean_pandas.py'),
                   '--input',
                   os.path.join(app.config['CLEAN_FOLDER'], code+'.gz'),
                   '--ntotal',
                   str(form.n_total.data),
                   '--ncas',
                   str(form.n_cases.data),
                   '--ncon',
                   str(form.n_controls.data)]
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            (output, err) = process.communicate()
            p_status = process.wait()
        else:
            flash('Error: size of output < 1MB')
    else:
        flash('Error: size of input < 1MB')

def munge_gwas(code):
    cmd = [os.path.join(app.config['PYTHON_LDSC_FOLDER'],'python'),
           os.path.join(app.config['LDSC_FOLDER'],'munge_sumstats.py'),
           '--out',
           os.path.join(app.config['MUNGED_FOLDER'], code),
           '--merge-alleles',
           os.path.join(app.config['LDSC_FOLDER'],'w_hm3.snplist'),
           '--N-cas-col',
           'Ncas',
           '--N-con-col',
           'Ncon',
           '--N-col',
           'N',
           '--sumstats',
           os.path.join(app.config['CLEAN_FOLDER'], code+'.gz'),
           '--info-min',
           '0.6']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    (output, err) = process.communicate()
    p_status = process.wait()                                         
    flash(output.decode())

def ldsc(code):
    cmd = [os.path.join(app.config['PYTHON_LDSC_FOLDER'],'python'),
           os.path.join(app.config['LDSC_FOLDER'],'ldsc.py'),
           '--out',
           os.path.join(app.config['MUNGED_FOLDER'], code+'_herit'),
           '--h2',
           os.path.join(app.config['MUNGED_FOLDER'], code+'.sumstats.gz'),
           '--ref-ld-chr',
           os.path.join(app.config['LDSC_FOLDER'],'eur_w_ld_chr/'),
           '--w-ld-chr',
           os.path.join(app.config['LDSC_FOLDER'],'eur_w_ld_chr/')]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    (output, err) = process.communicate()
    p_status = process.wait()                                         
    flash(output.decode())


