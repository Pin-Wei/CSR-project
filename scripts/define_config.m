function cfg = define_config()

cfg.do_norm = false; % true
cfg.norm_range = [-1, 1];

if cfg.do_norm
    cfg.out_tags = " (normed)";
else
    cfg.out_tags = "";
end

cfg.exclude_outlier = true; % false

if cfg.exclude_outlier
    cfg.out_tags = cfg.out_tags + " (out-rm)";
end

data_domains = ["psycholinguistic", "SRT"];
data_domain = data_domains(1);

switch data_domain

    case "psycholinguistic" % =============================================

        authors = ["Chang", "Tse"]; 
        author = authors(1); 
        
        tasks = ["Naming", "LD"]; 
        task = tasks(1); 

        cfg.data_folder = fullfile("..", "data", data_domain, join([author, task], "_")); 
    
        data_vers = ["raw", "zscored"];
        data_ver = data_vers(2);

        if data_ver == "zscored"
            cfg.out_tags = " (z-scored)" + cfg.out_tags;
        end
        
        data_levels = ["individual", "group", "single-subj"];
        data_level = data_levels(1);

        cfg.out_folder = fullfile("..", "output", data_domain, join([author, task], "_"), data_level); 
        
        switch data_level 

            case "individual" % -------------------------------------------
                
                cfg.fn_regex = "sub_*.xlsx";
                cfg.sid_regex = "sub_%d.xlsx";

                if data_ver == "zscored"
                    cfg.fn_regex = join(["zscored", cfg.fn_regex], "_");
                    cfg.sid_regex = join(["zscored", cfg.sid_regex], "_");
                end

            case "group" % ------------------------------------------------

                if data_ver == "raw"
                    cfg.fn_regex = "all_subjs_raw_*.xlsx";
                    cfg.sid_regex = "all_subjs_raw_%d (*).xlsx";
                else
                    cfg.fn_regex = "all_subjs_zvars_*.xlsx";
                    cfg.sid_regex = "all_subjs_zvars_%d (*).xlsx";
                end

            case "single-subj" % ------------------------------------------

                cfg.data_folder = fullfile(cfg.data_folder, "derivatives");
                sidx = 5102; % input("Which participant? ", "s");
                seed = 6180; % input("What was the seed? ", "s");

                cfg.fn_regex = join(["sub-", sidx, "_*.xlsx"], "_");
                cfg.sid_regex = join(["sub-", sidx, "_seed=", seed, "_%d.xlsx"], "_");

                if data_ver == "zscored"
                    cfg.fn_regex = join(["zscored", cfg.fn_regex], "_");
                    cfg.sid_regex = join(["zscored", cfg.sid_regex], "_");
                end

                cfg.out_tags = cfg.out_tags + " (sub-" + sidx + "_" + seed + ")";
        end

    case "SRT" % ==========================================================

        authors = ["BinShan", "Joanne"]; 
        author = authors(1); 

        switch author 

            case "BinShan" % ----------------------------------------------

                exp_names = ["Exp1_Spatial", "Exp2_Temporal", "Exp3_Intertwine(tDCS)", "Exp4_Intertwine"]; 
                exp_name = exp_names(1);
                
                data_vers = [".", "Basic_Shannon (2-back)", "Basic_Shannon (3-back)"];
                data_ver = data_vers(3);

                cfg.data_folder = fullfile("..", "data", data_domain, author, exp_name, "PreprocessedData", data_ver); 
                cfg.out_folder = fullfile("..", "output", data_domain, author, exp_name, data_ver); 

                if data_ver == "."
                    cfg.fn_regex = '*_Spatial_HL.csv';
                    cfg.sid_regex = '%d_Spatial_HL.csv';
                else    
                    cfg.fn_regex = 'sub_*.xlsx';
                    cfg.sid_regex = 'sub_%d.xlsx';
                end

            case "Joanne" % -----------------------------------------------

                data_vers = ["Basic_Shannon (2-back)", "Basic_Shannon (3-back)"];
                data_ver = data_vers(2);

                inst_types = ["explicit", "implicit"]; 
                inst_type = inst_types(1);
                
                subj_groups = ["learned", "not_learned"];
                subj_group = subj_groups(1);
                
                cfg.data_folder = fullfile("..", "data", data_domain, author, "PreprocessedData", data_ver, inst_type, subj_group); 
                cfg.out_folder = fullfile("..", "data", data_domain, author, data_ver, join([inst_type, subj_group], "_"));
                
                cfg.fn_regex = 'sub_*.xlsx';
                cfg.sid_regex = 'sub_%d.xlsx';
        end
end