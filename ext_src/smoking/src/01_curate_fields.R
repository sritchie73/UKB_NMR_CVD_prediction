library(data.table)

system("mkdir -p data/raw/ukbiobank/extracted/", wait=TRUE)

info <- rbind(use.names=FALSE,
	data.table(field.id=20160, name="Ever smoked", var="ever_smoked"),
	data.table(field.id=20162, name="Pack years adult smoking as proportion of life span exposed to smoking", var="pack_years_pct_lifespan"),
	data.table(field.id=20161, name="Pack years of smoking", var="pack_years"),
	data.table(field.id=10895, name="Light smokers, at least 100 smokes in lifetime (pilot)", var="light_smoker_pilot"),
	data.table(field.id=20116, name="Smoking status", "smoking_status"),
	data.table(field.id=1239, name="Current tobacco smoking", "current_tobacco"),
	data.table(field.id=1249, name="Past tobacco smoking", "past_tobacco"),
	data.table(field.id=2644, name="Light smokers, at least 100 smokes in lifetime", "light_smoker"),
	data.table(field.id=3436, name="Age started smoking in current smokers", "age_smoking_start"),
	data.table(field.id=3446, name="Type of tobacco currently smoked", "type_smoking"),
	data.table(field.id=5959, name="Previously smoked cigarettes on most/all days", "daily_past_smoker"),
	data.table(field.id=3456, name="Number of cigarettes currently smoked daily (current cigarette smokers)", "daily_cigarettes"),
	data.table(field.id=6194, name="Age stopped smoking cigarettes (current cigar/pipe or previous cigarette smoker)", "age_switch_to_cigar"),
	data.table(field.id=6183, name="Number of cigarettes previously smoked daily (current cigar/pipe smokers)", "past_daily_cigarettes_now_cigar"),
	data.table(field.id=3466, name="Time from waking to first cigarette", "wake_to_cigarette"),
	data.table(field.id=3476, name="Difficulty not smoking for 1 day", "difficulty_skip_day"),
	data.table(field.id=3486, name="Ever tried to stop smoking", "tried_quitting"),
	data.table(field.id=3496, name="Wants to stop smoking", "wants_to_quit"),
	data.table(field.id=3506, name="Smoking compared to 10 years previous", "current_vs_past10years"),
	data.table(field.id=6158, name="Why reduced smoking", "why_reduced_smoking"),
	data.table(field.id=2867, name="Age started smoking in former smokers", "age_start_after_stop"),
	data.table(field.id=2877, name="Type of tobacco previously smoked", "previous_type_smoking"),
	data.table(field.id=2887, name="Number of cigarettes previously smoked daily", "past_daily_cigarettes"),
	data.table(field.id=2897, name="Age stopped smoking", "quit_age"),
	data.table(field.id=2907, name="Ever stopped smoking for 6+ months", "ever_quit_6mo"),
	data.table(field.id=10827, name="Ever stopped smoking for 6+ months (pilot)", "ever_quit_6mo_pilot"),
	data.table(field.id=6157, name="Why stopped smoking", "why_stopped_smoking"),
	data.table(field.id=10115, name="Why stopped smoking (pilot)", "why_stopped_smoking_pilot"),
	data.table(field.id=2926, name="Number of unsuccessful stop-smoking attempts", "number_quit_attempts"),
	data.table(field.id=2936, name="Likelihood of resuming smoking", "likelihood_resume"),
	data.table(field.id=1259, name="Smoking/smokers in household", "household_smokers"),
	data.table(field.id=1269, name="Exposure to tobacco smoke at home", "home_smoke_exposure"),
	data.table(field.id=1279, name="Exposure to tobacco smoke outside home", "out_smoke_exposure")
)

fwrite(info, sep="\t", quote=FALSE, file="data/raw/ukbiobank/extracted/field_info.txt")
fwrite(info[,.(field.id)], quote=FALSE, col.names=FALSE, file="data/raw/ukbiobank/extracted/fields.txt")
