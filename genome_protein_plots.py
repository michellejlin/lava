#!/usr/bin/env pythonw

import argparse
import sys
import numpy as np
import pandas as pd
from bokeh.io import export_png, save
from bokeh.resources import CDN
from bokeh.embed import components, file_html, autoload_static
from bokeh.plotting import figure, show, output_file
from bokeh.palettes import Dark2, Category20, Category20c, Spectral10, Colorblind
from bokeh.models import ColumnDataSource, Jitter, BoxAnnotation, DataRange1d, HoverTool, CategoricalTicker, Slider, CustomJS, Label, WheelZoomTool, ResetTool, Button, TextInput
from bokeh.models.tickers import FixedTicker
from bokeh.transform import jitter, factor_cmap, linear_cmap
from bokeh.models.widgets import Panel, Tabs, Paragraph, Div, CheckboxGroup
from bokeh.layouts import column, layout, widgetbox, row
from palette import color_palette
import subprocess


#Uses info from protein csv to shade every other protein in light green, and annotate with protein names
def protein_annotation(first):
	protein_locs = []
	protein_names = []	
	#shades in every other protein region

	for i in range(0, proteins.shape[0]):
		if(i==0):
			x1 = 0
		elif(proteins.iloc[i,1] < proteins.iloc[i-1,2]) and first:
			print('WARNING: Protein-' + str(proteins.iloc[i,0]) + ' is overlapping with Protein-' + str(proteins.iloc[i-1,0]))
			print('Analysis will continue but the visualization for these two proteins will look a little funny. Often the fix for this is simply deleting the small ancilary proteins that overlapping from the gff file and using the -f and -g flags. For more help see the readme.')
			x1 = proteins.iloc[i,1]
		else:
			x1 = proteins.iloc[i,1]
		if(i==proteins.shape[0]-1):
			x2 = proteins.iloc[i,2]
		else:
			x2 = proteins.iloc[(i+1),1]
		if(i%2==0):
			g.add_layout(BoxAnnotation(left=x1, right=x2, fill_alpha=0.1, fill_color='green'))
		protein_locs.append((x1+x2)/2)
		protein_names.append(proteins.iloc[i,0])


	#adds protein labels as tick marks
	g.xaxis.ticker = protein_locs
	g.xaxis.major_label_overrides = dict(zip(protein_locs, protein_names))


#Creates the legend and configures some of the toolbar stuff
def configurePlot():
	g.legend.location = "top_right"
	g.legend.glyph_width = 35
	g.legend.glyph_height = 35
	g.legend.spacing = -10
	g.legend.background_fill_alpha = 0.5
	# disabled scroll wheel zooming for all plots - RCS
	#g.toolbar.active_scroll = g.select_one(WheelZoomTool)
	g.yaxis.axis_label_standoff = 10
	g.yaxis.axis_label = "Allele Frequency (%)"
	g.xaxis.minor_tick_line_color = None

#Callback for the sliders and the checkboxes for non/synonymous mutations, goes through data and selects/pushes the data with the matching criteria
def sliderCallback(ref_source, source, slider, slider_af, syngroup):
	return CustomJS(args=dict(ref=ref_source, depth_sample=source, slider=slider, slider_af=slider_af, checkbox=syngroup), code = """
		let depth = slider.value;
		let af = slider_af.value;
		let colnames = Object.keys(ref.data);
		for (let a = 0; a < colnames.length; a++) {
			depth_sample.data[colnames[a]] = [];
		}
		if (checkbox.active.indexOf(0)>-1) {
			for (let i = 0; i < ref.data['Depth'].length; i++) {
				if (depth <= ref.data['Depth'][i] && af <= ref.data['AF'][i] && ref.data['Syn'][i]=='synonymous SNV') {
					for (let b = 0; b < colnames.length; b++) {
						depth_sample.data[colnames[b]].push(ref.data[colnames[b]][i]);
					}
				}
			}
		}
		if (checkbox.active.indexOf(1)>-1) {
			for (let i = 0; i < ref.data['Depth'].length; i++) {
				if (depth <= ref.data['Depth'][i] && af <= ref.data['AF'][i] && ref.data['Syn'][i]=='nonsynonymous SNV') {
					for (let b = 0; b < colnames.length; b++) {
						depth_sample.data[colnames[b]].push(ref.data[colnames[b]][i]);
					}
				}
			}
		}
		if (checkbox.active.indexOf(2)>-1) {
			for (let i = 0; i < ref.data['Depth'].length; i++) {
				if (depth <= ref.data['Depth'][i] && af <= ref.data['AF'][i] && (ref.data['Syn'][i]=='stopgain'||ref.data['Syn'][i]=='stoploss')) {
					for (let b = 0; b < colnames.length; b++) {
						depth_sample.data[colnames[b]].push(ref.data[colnames[b]][i]);
					}
				}
			}
		}
		if (checkbox.active.indexOf(3)>-1) {
			for (let i = 0; i < ref.data['Depth'].length; i++) {
				if (depth <= ref.data['Depth'][i] && af <= ref.data['AF'][i] && (ref.data['Syn'][i]=='complex')) {
					for (let b = 0; b < colnames.length; b++) {
						depth_sample.data[colnames[b]].push(ref.data[colnames[b]][i]);
					}
				}
			}
		}
        depth_sample.change.emit();
    """)

# Second slider callback so that manual depth input overrides the slider 
def sliderCallback2(ref_source, source, slider, slider_af, syngroup, ose):
		return CustomJS(args=dict(ref=ref_source, depth_sample=source, slider=slider, slider_af=slider_af, checkbox=syngroup, ose=ose), code = """
			let depth = Number(ose.value);
			let af = slider_af.value;
			let colnames = Object.keys(ref.data);
			for (let a = 0; a < colnames.length; a++) {
				depth_sample.data[colnames[a]] = [];
			}
			if (checkbox.active.indexOf(0)>-1) {
				for (let i = 0; i < ref.data['Depth'].length; i++) {
					if (depth <= ref.data['Depth'][i] && af <= ref.data['AF'][i] && ref.data['Syn'][i]=='synonymous SNV') {
						for (let b = 0; b < colnames.length; b++) {
							depth_sample.data[colnames[b]].push(ref.data[colnames[b]][i]);
						}
					}
				}
			}
			if (checkbox.active.indexOf(1)>-1) {
				for (let i = 0; i < ref.data['Depth'].length; i++) {
					if (depth <= ref.data['Depth'][i] && af <= ref.data['AF'][i] && ref.data['Syn'][i]=='nonsynonymous SNV') {
						for (let b = 0; b < colnames.length; b++) {
							depth_sample.data[colnames[b]].push(ref.data[colnames[b]][i]);
						}
					}
				}
			}
			if (checkbox.active.indexOf(2)>-1) {
				for (let i = 0; i < ref.data['Depth'].length; i++) {
					if (depth <= ref.data['Depth'][i] && af <= ref.data['AF'][i] && (ref.data['Syn'][i]=='stopgain'||ref.data['Syn'][i]=='stoploss')) {
						for (let b = 0; b < colnames.length; b++) {
							depth_sample.data[colnames[b]].push(ref.data[colnames[b]][i]);
						}
					}
				}
			}
			if (checkbox.active.indexOf(3)>-1) {
				for (let i = 0; i < ref.data['Depth'].length; i++) {
					if (depth <= ref.data['Depth'][i] && af <= ref.data['AF'][i] && (ref.data['Syn'][i]=='complex')) {
						for (let b = 0; b < colnames.length; b++) {
							depth_sample.data[colnames[b]].push(ref.data[colnames[b]][i]);
						}
					}
				}
			}
            depth_sample.change.emit();
        """)
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='no description yet')
	parser.add_argument('-png', action='store_true')
	parser.add_argument('-nuc', action='store_true')
	parser.add_argument('-pdb', help='Provide a PDB acsession number to load into the NGL viewer')
	parser.add_argument('merged', help='internal call pointing to the merged.csv file created by the main lava script')
	parser.add_argument('proteins', help='Internal call pointing to the proteins.csv file created by the main lava script')
	parser.add_argument('reads', help='Internal call pointing to the reads.csv file created by the main lava script')
	parser.add_argument('dir', help='Dir created by main lava program')
	
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)
	pdb_num = args.pdb 
	new_dir = args.dir

	# read data from lava output
	merged = pd.read_csv(args.merged,index_col=False)
	proteins = pd.read_csv(args.proteins,index_col=False,header=None)
	reads = pd.read_csv(args.reads, index_col=False, names = ["Sample", "Total", "Mapped", "Percentage"])

	source = ColumnDataSource(merged)
	list_tabs = []
	list2_tabs = []
	list_plots = []
	list2_plots = []

	num_samples = merged.Sample.nunique()
	unique_samples = merged.Sample.unique()
	unique_samples = np.sort(merged.Sample.unique())
	num_Passages = merged['Passage'].max()

	TOOLTIPS = [
		("Amino Acid Change", "@Change"),
		("Nucleotide Change", "@NucleotideChange"),
		("Allele Frequency", "@AF"+"%"),
		("Depth", "@Depth"),
	]

	#creates genome plots
	FIRST = True
	# go through each sample and create the genome plot 
	for sample_name, merged_Sample in merged.groupby('Sample', sort=False):
		sample = merged_Sample
		name = sample_name
		
		# Plot per base coverage for each sample as a subplot 
		coverage = pd.read_table(sample_name.strip() + '.genomecov', names=["sample", 'position', 'cov'], header=0)
		cov_sample=ColumnDataSource(coverage)
		cov_sample.data['position'].astype(float)
		cov_val_sample = ColumnDataSource(data=cov_sample.data)
		
		f = figure(plot_width=400, plot_height=200, title='Coverage')
		f.line(x='position', y='cov',source=cov_val_sample)
		
		# setup 
		source_sample=ColumnDataSource(merged_Sample)
		source_sample.data['Position'].astype(float)
		depth_sample = ColumnDataSource(data=source_sample.data)
		#creates a graph with xaxis being genome length (based on protein csv)
		# took out active_scroll = "wheel_zoom" -RCS
		g = figure(plot_width=1600, plot_height=800, y_range=DataRange1d(bounds=(0,102), start=0,end=102),
			title=sample_name.split('/')[1], sizing_mode = 'scale_width',
			x_range=DataRange1d(bounds=(0, proteins.iloc[proteins.shape[0]-1,2]), start=0, end=proteins.iloc[proteins.shape[0]-1,2]))
		#graphs scatterplot, with different colors for non/synonymous mutations
		# this creates console errors when complex mutations are detected, can probobly fix this by adding more colors 
		if(args.nuc):
			g.circle(x='Position', y=jitter('AF', width=2, range=g.y_range), size=15, alpha=0.6, hover_alpha=1, 
				legend = 'LetterChange', line_color='white', line_width=2, line_alpha=1,
				fill_color=factor_cmap('LetterChange', palette=color_palette, factors=merged.LetterChange.unique()), 
				hover_color=factor_cmap('LetterChange', palette=color_palette, factors=merged.LetterChange.unique()),
				source=depth_sample, hover_line_color='white')
		else:
			g.circle(x='Position', y=jitter('AF', width=2, range=g.y_range), size=15, alpha=0.6, hover_alpha=1, 
				legend = 'Syn', line_color='white', line_width=2, line_alpha=1,
				fill_color=factor_cmap('Syn', palette=['firebrick', 'cornflowerblue', 'green', 'purple', 'yellow'], factors=merged.Syn.unique()), 
				hover_color=factor_cmap('Syn', palette=['firebrick', 'cornflowerblue', 'green', 'purple', 'yellow'], factors=merged.Syn.unique()),
				source=depth_sample, hover_line_color='white')
		
		g.add_tools(HoverTool(tooltips=TOOLTIPS))
		g.xgrid.grid_line_alpha = 0
		configurePlot()
		protein_annotation(FIRST)
		FIRST = False
		g.xaxis.axis_label = "Protein"

		b = Button(label="Reset Plot")
		b.js_on_click(CustomJS(args=dict(g=g), code="""
			g.reset.emit()
		"""))
		ose = TextInput(title='Manually input depth:')
		
		syngroup = CheckboxGroup(labels=["Show synonymous mutations", "Show nonsynonymous mutations", 
			"Show stopgains and stoplosses", "Show complex mutations"], active=[0,1,2,3])
		
		slider = Slider(start=0, end=merged['Depth'].max(), step=1, value=5, title='Depth Cutoff')
		slider_af = Slider(start=0, end=100, step=1, value=5, title='Allele Frequency Cutoff')
		slider.js_on_change('value', sliderCallback(source_sample, depth_sample, slider, slider_af, syngroup))
		slider_af.js_on_change('value', sliderCallback(source_sample, depth_sample, slider, slider_af, syngroup))
		syngroup.js_on_change('active', sliderCallback(source_sample, depth_sample, slider, slider_af, syngroup))
		ose.js_on_change('value', sliderCallback2(source_sample, depth_sample, slider, slider_af, syngroup, ose ))
		
		#labels with read information
		reads_info = (reads.loc[reads['Sample'] == name.strip()])
		div=Div(text="""<font size="2" color="gray"><b>Total reads: </b>"""+str(reads_info.iloc[0]["Total"])+
			"""<br><b>Total reads mapped: </b>"""+str(reads_info.iloc[0]["Mapped"])+"""<br><b>Percentage of reads mapped: </b>
			"""+str(reads_info.iloc[0]["Percentage"])+"""</font>""", width=300, height=100)
		
		g = layout(row([g, column([Div(text="""""", width=300, height=220),f,ose, slider, slider_af, syngroup, widgetbox(div),b])]))
    	
		#creates both tabs and different plots
		tab = Panel(child=g, title=name.split('/')[1])
		list_tabs.append(tab)
		list_plots.append(g)
	
	tabs_genomes = Tabs(tabs=list_tabs)
	plots_genomes = (column(list_plots))
	
	total_num_mutations = 0
	#creates protein plots
	for protein_name, merged_Protein in merged.groupby('Protein'):
		protein = merged_Protein
		name = protein_name
		
		source_protein = ColumnDataSource(protein)
		depth_sample_p = ColumnDataSource(data=source_protein.data)
		#CHANGE RANGE AS NEEDED
		g = figure(plot_width=1600, plot_height=800, y_range=DataRange1d(bounds=(0,102), start=0,end=102),
			title=protein_name,
 			x_range=DataRange1d(bounds=(-1,num_Passages+2), start=-1, end=num_Passages+2))
			#x_range=DataRange1d(bounds=(-1,5), start=-1, end=5))
		
		#creates list of lists by getting all the xvalues and yvalues for multiline
		xlist_all = []
		ylist_all = []
		protein_mut = merged.loc[merged['Protein']==protein_name]
		
		num_mutations = 0
		unique_mutations = (protein_mut.groupby('Change')[['Passage', 'AF']].apply(lambda x: x.values.tolist()))
		for unique_mutation in unique_mutations:
			xlist = []
			ylist = []
			for individual in unique_mutation:
				xlist.append(individual[0])
				ylist.append(individual[1])
			xlist_all.append(xlist)
			ylist_all.append(ylist)
			num_mutations+=1
		
		#makes visible multiline when hovered over
		#g.multi_line(xlist_all, ylist_all, color=palette[total_num_mutations:(total_num_mutations+num_mutations)])
		line = g.multi_line(xlist_all, ylist_all, line_width=2, alpha=0, hover_line_alpha=0.6,
			hover_line_color="gray", line_cap = 'round', line_dash = 'dotted')
		
		total_num_mutations+=num_mutations
		
#		circle = g.circle(x=jitter('Passage', width=0.3, range=g.x_range), y=jitter('AF', width=2, range=g.y_range), size=12, alpha=0.8, 
		if (args.png):
			circle = g.circle(x=jitter('Passage',width=0.35, range=g.x_range), y='AF', size=15, alpha=0.8,
				fill_color=factor_cmap('Change', palette=color_palette, factors=merged.Change.unique()),
				line_color='white', line_width=2, line_alpha=1,legend = 'Change',source=depth_sample_p)
		else:
			circle = g.circle(x=jitter('Passage',width=0.08), y='AF', size=15, alpha=0.7, hover_alpha = 1,
				fill_color=factor_cmap('Change', palette=color_palette, factors=merged.Change.unique()), 
				hover_color=factor_cmap('Change', palette=color_palette, factors=merged.Change.unique()),
				line_color='white', line_width=2, hover_line_color='white', line_alpha=1,legend = 'Change',source=depth_sample_p)

		g.xaxis.axis_label = "Passage"
		#g.xaxis.ticker = FixedTicker(ticks=[0,1,3,5])
		
		configurePlot()
		#don't want the tooltips to show up for multiline, so makes separate hovertools for each glyph
		ht = HoverTool(renderers=[circle], tooltips=TOOLTIPS)
		ht2 = HoverTool(renderers=[line], tooltips=None)
		g.tools.append(ht)
		g.tools.append(ht2)
		
		b = Button(label="Reset Plot")
		b.js_on_click(CustomJS(args=dict(g=g), code="""
			g.reset.emit()
		"""))
		
		ose = TextInput(title='Manually input depth:')
		
		syngroup = CheckboxGroup(labels=["Show synonymous mutations", "Show nonsynonymous mutations", 
			"Show stopgains and stoplosses", "Show complex mutations"], active=[0,1,2,3])
					
		slider = Slider(start=0, end=merged['Depth'].max(), step=1, value=5, title='Depth Cutoff')
		slider_af = Slider(start=0, end=100, step=1, value=5, title='Allele Frequency Cutoff')
		slider.js_on_change('value', sliderCallback(source_protein, depth_sample_p, slider, slider_af,syngroup))
		slider_af.js_on_change('value', sliderCallback(source_protein, depth_sample_p, slider, slider_af,syngroup))
		syngroup.js_on_change('active', sliderCallback(source_protein, depth_sample_p, slider, slider_af, syngroup))
		ose.js_on_change('value', sliderCallback2(source_protein, depth_sample_p, slider, slider_af, syngroup, ose ))
		g = layout(row([g, column([Div(width=100, height=260),ose,slider, slider_af, syngroup, Div(width=100, height=30),b])]))
		tab = Panel(child=g, title=name)
		list2_tabs.append(tab)
		list2_plots.append(g)		
	
	tabs_proteins = Tabs(tabs=list2_tabs)
	plots_proteins = (column(list2_plots))
	
#can only either have html or pngs because of weird Bokeh voodoo
	if(args.nuc):
		output_file("nucleotide_changes.html")
		show(tabs_genomes)
	else:
		if(args.png):
			export_png(plots_proteins, filename="Protein_Plots.png")
			export_png(plots_genomes, filename="Genome_Plots.png")
		else:
			# save output both as standalone HTML file and as a javascript element and a script tag 
			print('Opening output file genome_protein_plots.html...\nGraphs_and_viewer.html includes the protein viewer.')
			output_file(new_dir + "/" + new_dir + "_plots.html")
			subprocess.call('cp ngls_test.html ' + new_dir + '/', shell=True)
			save(column(tabs_genomes, tabs_proteins))
			# now the file exists
			plot_file = open(new_dir + '/' + new_dir + '_plots.html')
			viewer_code_file = open(new_dir + '/' + 'ngls_test.html')
			new_file = open(new_dir + '/' + new_dir + '_plots_and_viewer.html', 'w')
			script_tag_count = 0
			for line in plot_file:
				new_file.write(line)
				if line.strip() == '</script>':
					script_tag_count += 1 
				if script_tag_count == 3:
					for line_two in viewer_code_file:
						new_file.write(line_two)
			plot_file.close()
			viewer_code_file.close()
			new_file.close()
			try:
				show(column(tabs_genomes, tabs_proteins))
			except:
				print('Automatic opening of output files has failed - however generation worked. Check out your output at the above paths.')
				

