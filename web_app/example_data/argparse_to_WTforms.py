#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: vbhardwaj
"""
import sys
import re
import importlib.util
import click

# all argparse keys (I think)
ALLKEYS={'choices',
 'const',
 'container',
 'default',
 'dest',
 'help',
 'metavar',
 'nargs',
 'option_strings',
 'required',
 'type'}

def arg_to_field(arg):
    argkeys=set(arg.keys())
    newkeys=argkeys - ALLKEYS
    if len(newkeys) > 0:
        print("Some keys in arg are not in expected keys: \n")
        print(newkeys)

    arg_dest=arg['dest']
    arg_name=re.sub('([A-Z][a-z]+)', r' \1', re.sub('([A-Z]+)', r' \1', re.sub("--", "", arg['option_strings'][0]))).title()
    arg_default=arg['default']
    arg_help=re.sub(' \(Default: \%\(default\)s\)', '', arg['help'])

    # identify the type of field
    if arg['choices']:
            arg_field='SelectField'
    elif arg['type']:
        if callable(arg['type']):
            if arg['type'] is int:
                arg_field='IntegerField'
            elif arg['type'] is float:
                arg_field='FloatField'
            elif arg['type'] is str:
                arg_field='StringField'
            else:
                arg_field='FileField'

    else:
        if arg['metavar'] == 'FILE':
            if arg['nargs'] == '+':
                arg_field='MultipleFileField'
            else:
                arg_field='FileField'
        else:
            arg_field='StringField'


    # add validators
    validator_list=[]
    if arg['required']:
        validator_list.append('validators.DataRequired()')
    arg_validator=",".join(validator_list)

    ## print a string based on the information here
    if arg_field=='SelectField':
        arg_choices=arg['choices']
        out = """{}={}(label='{}', validators=[{}],
                 choices='{}', default='{}', description='{}')""".format(arg_dest,
                                                                            arg_field,
                                                                            arg_name,
                                                                            arg_validator,
                                                                            arg_choices,
                                                                            arg_default,
                                                                            arg_help)
    else:
        out = "{}={}(label='{}', validators=[{}], default='{}', description='{}')".format(arg_dest,
                                                                            arg_field,
                                                                            arg_name,
                                                                            arg_validator,
                                                                            arg_default,
                                                                            arg_help)
    return out


@click.command()
@click.option('--inargs', '-a')
@click.option('--infile', '-f')
def run(inargs, infile):

    spec = importlib.util.spec_from_file_location("module.name", infile)
    foo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(foo)
    p=foo.parseArguments()
    p_out=p.parse_args(inargs.split())

    ## print to stdout
    for ac in p.__dict__['_actions']:
        info=vars(ac)
        if info['dest']=='help' or info['dest']=='verbose' or info['dest']=='version':
            continue
        print(arg_to_field(info))
    print("submit = SubmitField('Submit')")

if __name__ == "__main__":
    run()
