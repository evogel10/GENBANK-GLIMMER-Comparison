#!/usr/local/bin/python3

import jinja2
import re

# This line tells the template loader where to search for template files
templateLoader = jinja2.FileSystemLoader( searchpath="./templates" )

# This creates your environment and loads a specific template
env = jinja2.Environment(loader=templateLoader)
template = env.get_template('comparison.html')

genbank_file = open('sequence.gb', 'r')
glimmer_file = open('run1.predict', 'r')

# Tracks Reference Gene/CDS entries
genbank_count = 0
# Tracks Glimmer entries
glimmer_count = 0
# Tracks exact matches
exact_count = 0
# Tracks 5 prime matches
five_count = 0
# Tracks 3 prime matches
three_count = 0
# Tracks no matches
no_count = 0
# Flag to see if a new CDS has been reached
new_cds = False
# Stores GenBank gene ID and start/stop positions
genbank_info = []
# Stores Glimmer gene ID and start/stop positions
glimmer_info = []
# GenBank start position
start = 0
# GenBank stop position
stop = 0
# GenBank ID
ref_id = ''
# GenBank entry
gen = -1
# Glimmer entry
glim = -1
#tracks data to push to data file
data_table = []


# Finds gene ID and start/stop positions from GenBank file
for line in genbank_file:
	# Removes whitespce on both sides of string
	line = line.strip()
	# Splits the string into a list on whitespaces
	line = line.split(' ')
	# Removes empty space inside of the list
	line = list(filter(None, line))
	# Holds gene ID and start/stop positions
	gen_temp_list = []

	# Finds the CDS regions
	if len(line) > 0:

		if line[0] == 'CDS':
			genbank_count += 1
			positions = line[1].split('..')

			# Look for the first postions with the join before it
			if line[1][0] == 'j':
				start = re.sub("\D", "", positions[0])
				stop = re.sub("\D", "", positions[2])
			# finds complement positions
			elif line[1][0] == 'c':
				start_temp = re.sub("\D", "", positions[1])
				start = '-' + start_temp
				stop_temp = re.sub("\D", "", positions[0])
				stop = '-' + stop_temp
			# finds the standard 5' 3' positions
			elif line[1][0].isdigit():
				start = positions[0]
				stop = positions[1]
			# finds the positions with < and > infront
			elif line[1][0] == '>' or line[1][0] == '<':
				start = positions[0][1:]
				stop = positions[1]

			# Marker to show which CDS you are getting information from
			new_cds = True
	
	# Finds the genbank_id
	if new_cds == True:	
		if line[0].startswith('/protein_id='):
			ids = line[0].split('"')
			ref_id = ids[1]
			new_cds = False
		if ref_id != '':
			gen_temp_list.append(ref_id)
			gen_temp_list.append(start)
			gen_temp_list.append(stop)
			genbank_info.append(gen_temp_list)
			ref_id = ''

# Finds gene ID and start/stop positions from Glimmer3 file
for line in glimmer_file:
	# Holds gene ID and start/stop positions
	glim_temp_list = []

	# Skips the first line in the file describing the sequence
	if line.startswith('>'):
		continue

	glimmer_count += 1

	glimmer_line = line.strip().split(' ')
	# Removes empty space inside of the list
	glimmer_line = list(filter(None, glimmer_line))

	# Finds glimmer_id
	glim_temp_list.append(glimmer_line[0])

	# adds the positions while finding complements
	if glimmer_line[3][0] == '+':
		glim_temp_list.append(glimmer_line[1])
		glim_temp_list.append(glimmer_line[2])
	elif glimmer_line[3][0] == '-':

		glim_temp_list.append('-' + glimmer_line[1])
		glim_temp_list.append('-' + glimmer_line[2])

	glimmer_info.append(glim_temp_list)

# Finds exact, 5 prime, 3 prime, and no matches between GenBank and Glimmer files
while gen < (len(genbank_info) - 1) and glim < (len(glimmer_info) - 1):
	temp_data_table = []

	# Iterates to next GenBank entry
	gen +=1
	# Iterates to next Glimmer entry
	glim +=1
	# 5 prime flag
	match_5 = False
	# 3 prime flag
	match_3 = False
	# Flag that shows if there was a match after iterating through Glimmer file
	match_flag = False

	# Looks for 5 prime match
	if genbank_info[gen][1] == glimmer_info[glim][1]:
		match_5 = True

	# Looks for 3 prime match
	if genbank_info[gen][2] == glimmer_info[glim][2]:
		match_3 = True

	# Checks for matches
	if match_5 == True and match_3 == True:

		# print(genbank_info[gen][0] + ' | ' + genbank_info[gen][1] + ' | ' + genbank_info[gen][2] + ' | ' + glimmer_info[glim][0] + ' | ' + glimmer_info[glim][1] + ' | ' + glimmer_info[glim][2] + ' | EXACT' )
		temp_data_table.append(genbank_info[gen][0])
		temp_data_table.append(genbank_info[gen][1])
		temp_data_table.append(genbank_info[gen][2])
		temp_data_table.append(glimmer_info[glim][0])
		temp_data_table.append(glimmer_info[glim][1])
		temp_data_table.append(glimmer_info[glim][2])
		temp_data_table.append('EXACT')
		exact_count += 1

	elif match_5 == True and match_3 == False:

		# print(genbank_info[gen][0] + ' | ' + genbank_info[gen][1] + ' | ' + genbank_info[gen][2] + ' | ' + glimmer_info[glim][0] + ' | ' + glimmer_info[glim][1] + ' | ' + glimmer_info[glim][2] + ' | 5 MATCH' )
		temp_data_table.append(genbank_info[gen][0])
		temp_data_table.append(genbank_info[gen][1])
		temp_data_table.append(genbank_info[gen][2])
		temp_data_table.append(glimmer_info[glim][0])
		temp_data_table.append(glimmer_info[glim][1])
		temp_data_table.append(glimmer_info[glim][2])
		temp_data_table.append('5 MATCH')
		count += 1
		five_count += 1

	elif match_5 == False and match_3 == True:

		# print(genbank_info[gen][0] + ' | ' + genbank_info[gen][1] + ' | ' + genbank_info[gen][2] + ' | ' + glimmer_info[glim][0] + ' | ' + glimmer_info[glim][1] + ' | ' + glimmer_info[glim][2] + ' | 3 MATCH' )
		temp_data_table.append(genbank_info[gen][0])
		temp_data_table.append(genbank_info[gen][1])
		temp_data_table.append(genbank_info[gen][2])
		temp_data_table.append(glimmer_info[glim][0])
		temp_data_table.append(glimmer_info[glim][1])
		temp_data_table.append(glimmer_info[glim][2])
		temp_data_table.append('3 MATCH')
		three_count += 1

	else:

		whileflag = False
		# Track iterating through the glimmer file
		loop = 0

		# If no matches iterate through glimmer file
		while (gen < len(genbank_info) and glim < len(glimmer_info)) and ((abs(int(genbank_info[gen][1])) > abs(int(glimmer_info[glim][1])) and match_flag == False) or (((int(genbank_info[gen][1]) > 0) != ((int(glimmer_info[glim][1])) > 0)) and loop == 0)):

			whileflag = True
			glim += 1
			loop+=1
			if genbank_info[gen][1] == glimmer_info[glim][1]:
				match_5 = True
			if genbank_info[gen][2] == glimmer_info[glim][2]:
				match_3 = True
			if match_5 == True and match_3 == True:

				# print(genbank_info[gen][0] + ' | ' + genbank_info[gen][1] + ' | ' + genbank_info[gen][2] + ' | ' + glimmer_info[glim][0] + ' | ' + glimmer_info[glim][1] + ' | ' + glimmer_info[glim][2] + ' | EXACT - While Loop' )
				temp_data_table.append(genbank_info[gen][0])
				temp_data_table.append(genbank_info[gen][1])
				temp_data_table.append(genbank_info[gen][2])
				temp_data_table.append(glimmer_info[glim][0])
				temp_data_table.append(glimmer_info[glim][1])
				temp_data_table.append(glimmer_info[glim][2])
				temp_data_table.append('EXACT')
				match_flag = True
				exact_count += 1

			elif match_5 == True and match_3 == False:

				# print(genbank_info[gen][0] + ' | ' + genbank_info[gen][1] + ' | ' + genbank_info[gen][2] + ' | ' + glimmer_info[glim][0] + ' | ' + glimmer_info[glim][1] + ' | ' + glimmer_info[glim][2] + ' | 5 MATCH - While loop' )
				temp_data_table.append(genbank_info[gen][0])
				temp_data_table.append(genbank_info[gen][1])
				temp_data_table.append(genbank_info[gen][2])
				temp_data_table.append(glimmer_info[glim][0])
				temp_data_table.append(glimmer_info[glim][1])
				temp_data_table.append(glimmer_info[glim][2])
				temp_data_table.append('5 MATCH')
				five_count += 1
				match_flag = True

			elif match_5 == False and match_3 == True:

				# print(genbank_info[gen][0] + ' | ' + genbank_info[gen][1] + ' | ' + genbank_info[gen][2] + ' | ' + glimmer_info[glim][0] + ' | ' + glimmer_info[glim][1] + ' | ' + glimmer_info[glim][2] + ' | 3 MATCH - While loop' )
				temp_data_table.append(genbank_info[gen][0])
				temp_data_table.append(genbank_info[gen][1])
				temp_data_table.append(genbank_info[gen][2])
				temp_data_table.append(glimmer_info[glim][0])
				temp_data_table.append(glimmer_info[glim][1])
				temp_data_table.append(glimmer_info[glim][2])
				temp_data_table.append('3 MATCH')
				three_count += 1
				match_flag = True

			elif match_5 == False and match_3 == False:

				match_flag = False

				if(loop==1):
				 	glim -= 2

			match_3=False
			match_5=False

		# Tracks no matches
		if(match_flag == False and whileflag == False ):
			# print(genbank_info[gen][0] + ' | ' + genbank_info[gen][1] + ' | ' + genbank_info[gen][2] + ' | ' + glimmer_info[glim][0] + ' | ' + glimmer_info[glim][1] + ' | ' + glimmer_info[glim][2] + ' | NO Match' )
			temp_data_table.append(genbank_info[gen][0])
			temp_data_table.append(genbank_info[gen][1])
			temp_data_table.append(genbank_info[gen][2])
			temp_data_table.append(glimmer_info[glim][0])
			temp_data_table.append(glimmer_info[glim][1])
			temp_data_table.append(glimmer_info[glim][2])
			temp_data_table.append('EXACT')
			no_count+=1
			glim-=1
			if(((int(genbank_info[gen][1]) > 0) != ((int(glimmer_info[glim][1])) > 0))):
				glim-=1
		if(match_flag == False and whileflag == True ):
			# print(genbank_info[gen][0] + ' | ' + genbank_info[gen][1] + ' | ' + genbank_info[gen][2] + ' | ' + glimmer_info[glim][0] + ' | ' + glimmer_info[glim][1] + ' | ' + glimmer_info[glim][2] + ' | NO Match' )
			temp_data_table.append(genbank_info[gen][0])
			temp_data_table.append(genbank_info[gen][1])
			temp_data_table.append(genbank_info[gen][2])
			temp_data_table.append(glimmer_info[glim][0])
			temp_data_table.append(glimmer_info[glim][1])
			temp_data_table.append(glimmer_info[glim][2])
			temp_data_table.append('EXACT')
			no_count+=1
			glim-=1
			if(((int(genbank_info[gen][1]) > 0) != ((int(glimmer_info[glim][1])) > 0))):
				glim-=1

	data_table.append(temp_data_table)

print("Content-Type: text/html\n\n")
print(template.render(gen_count=genbank_count, glim_count=glimmer_count, exact=exact_count, five=five_count, three=three_count, no=no_count, data=data_table))



























