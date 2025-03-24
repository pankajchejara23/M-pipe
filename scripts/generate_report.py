import argparse
import os
from jinja2 import Environment, FileSystemLoader

def parse_arguments():
    """Parse command-line arguments for file paths."""
    parser = argparse.ArgumentParser(description="Generate an HTML microbiome report with given plot files.")
    parser.add_argument("--alpha_plot", required=True, help="Path to alpha diversity plot")
    parser.add_argument("--top_taxa_plot", required=True, help="Paths to top taxa plots (any number of files)")
    parser.add_argument("--diff_deseq2_plots", nargs='+', required=True, help="Path to differential analysis plot")
    #parser.add_argument("--differential_results", nargs='+', required=True, help="Path to differential analysis results file (CSV)")
    parser.add_argument("--auc_plot", required=True, help="Path to AUC plot for ML model performance")
    parser.add_argument("--auc_report", required=True, help="Path to model report ")
    parser.add_argument("--output", default="microbiome_report.html", help="Output HTML report file name")
    
    return parser.parse_args()

def generate_report(args):
    """Generate the HTML report using Jinja2 and user-provided file paths."""
    
    # Get the directory of the current script (scripts/)
    script_dir = os.path.dirname(__file__)

    # Point to the templates folder, which is above the script directory
    template_dir = os.path.join(script_dir, '..', 'template')


    # Load Jinja2 template
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template("base-template.html")

 

    # Define file paths from arguments
    plot_paths = {
        "alpha_plot": args.alpha_plot,
        "top_taxa_plot": args.top_taxa_plot,  
        "auc_plot": args.auc_plot,
        "auc_measure":.78,
        "diff_deseq2_plots":args.diff_deseq2_plots
        
    }

    # Render the template with actual plot paths
    rendered_html = template.render(plot_paths)

    # Save the final report
    with open(args.output, "w", encoding="utf-8") as f:
        f.write(rendered_html)

    print(f"Report generated: {args.output}")

if __name__ == "__main__":
    args = parse_arguments()
    print(args)
    generate_report(args)
