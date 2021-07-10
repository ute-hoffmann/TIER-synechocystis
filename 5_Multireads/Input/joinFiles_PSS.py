
first_file = "PSS_5ends_coverage_multireads/dIF1-1_0h_PSS.tabular"
files = ["dIF1-1_1h_PSS.tabular", "dIF1-2_0h_PSS.tabular", "dIF1-2_1h_PSS.tabular", "dIF3-1_0h_PSS.tabular", "dIF3-1_1h_PSS.tabular",  "dWT1-1_0h_PSS.tabular", "dWT1-1_1h_PSS.tabular", "dWT2_0h_PSS.tabular", "dWT2_1h_PSS.tabular", "dWT3_0h_PSS.tabular", "dWT3_1h_PSS.tabular"]
dir_files = []
direct = "PSS_5ends_coverage_multireads"
for i in files:
	dir_files.append(direct + "/" + i)

data_dict = {}
nt_list = []

with open(first_file) as dIF11_0h:
	print("now reading first file")
	for line in dIF11_0h:
		split_line = line.strip("\n").split("\t")
		nt_pos = split_line[0]
		count = split_line[1]
		nt_list.append(nt_pos)
		data_dict[nt_pos] = [count]
			
for tabular in dir_files:
	with open(tabular) as data:
		print("now reading " + tabular)
		for line in data:
			split_line = line.strip("\n").split("\t")
			nt_pos = split_line[0]
			count = split_line[1]
			if nt_pos not in data_dict.keys():
				print(nt_pos + " not in data dict !!!!")
			data_dict[nt_pos].append(count)
				
combined_table = "multireads_PSS_5ends_combined.txt"
with open(combined_table, "w") as table:
	print("now printing file")
	table.write("\tdIF11_0h\tdIF11_1h\tdIF12_0h\tdIF12_1h\tdIF31_0h\tdIF31_1h\tdWT1_0h\tdWT1_1h\tdWT2_0h\tdWT2_1h\tdWT3_0h\tdWT3_1h\n")
	for nt_pos in nt_list:
		table.write(nt_pos + "\t")
		j = 1
		for i in data_dict[nt_pos]:
			if j == 12:
				table.write(i + "\n")
				continue
			table.write(i + "\t")
			j += 1
