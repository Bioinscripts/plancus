import sys

def cluster_bed_to_gff(bed_sorted: str, gff_out: str, max_gap: int = 1000):
    clusters = []
    current_cluster = []
    
    with open(bed_sorted, 'r') as f:
        for line in f:
            if not line.strip(): continue
            parts = line.strip().split('\t')
            contig = parts[0]
            start = int(parts[1]) + 1 # BED (0-based) to GFF (1-based)
            end = int(parts[2])
            orf_id = parts[3]
            strand = parts[5] if len(parts) >= 6 else '+'
            
            entry = {"contig": contig, "start": start, "end": end, "strand": strand, "id": orf_id}
            
            if not current_cluster:
                current_cluster = [entry]
            else:
                last = current_cluster[-1]
                if contig == last["contig"] and (start - last["end"]) <= max_gap:
                    current_cluster.append(entry)
                else:
                    clusters.append(current_cluster)
                    current_cluster = [entry]
                    
        if current_cluster:
            clusters.append(current_cluster)
            
    bed_clusters_out = gff_out.replace(".gff3", "_clean.bed")
    with open(gff_out, 'w') as f, open(bed_clusters_out, 'w') as fb:
        f.write("##gff-version 3\n")
        cluster_id = 1
        for cl in clusters:
            c_contig = cl[0]["contig"]
            c_start = min(x["start"] for x in cl)
            c_end = max(x["end"] for x in cl)
            strands = [x["strand"] for x in cl]
            c_strand = max(set(strands), key=strands.count)
            
            cluster_name = f"{c_contig}_{c_start}_{c_end}_EVE_Cluster_{cluster_id}"
            
            f.write(f"{c_contig}\tEVE_Pipeline\tregion\t{c_start}\t{c_end}\t.\t{c_strand}\t.\tID={cluster_name};Name=Cluster_{cluster_id}\n")
            fb.write(f"{c_contig}\t{c_start-1}\t{c_end}\tCluster_{cluster_id}\t1000\t{c_strand}\n")
            
            orf_idx = 1
            for orf in cl:
                f.write(f"{c_contig}\tEVE_Pipeline\tCDS\t{orf['start']}\t{orf['end']}\t.\t{orf['strand']}\t.\tID={orf['id']};Parent={cluster_name};Name=ORF_{orf_idx}\n")
                orf_idx += 1
                
            cluster_id += 1

if __name__ == "__main__":
    if len(sys.argv) >= 4:
        cluster_bed_to_gff(sys.argv[1], sys.argv[2], int(sys.argv[3]))
