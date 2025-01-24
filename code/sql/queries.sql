-- Example of how ONS does January rebasing
SELECT date, cd.item_id, item_desc, item_index
FROM cpi_data cd
    INNER JOIN items i ON cd.item_id = i.item_id
WHERE cd.item_id = 610310
    AND EXTRACT(MONTH FROM date) IN (12, 1, 2)
    AND EXTRACT(YEAR FROM date) < 2011
ORDER BY date
OFFSET 1;

/*
Query to get the relevant data for analysis. Rebases the indices
so that they are relative to one base year.
*/
WITH first_years AS (
    -- Get the first year of data for each item
    SELECT
        item_id,
        MIN(EXTRACT(YEAR FROM date)) as first_year
    FROM cpi_data
    GROUP BY item_id
),
january_indices AS (
    -- Extract January indices and first year records
    SELECT
        cd.item_id,
        EXTRACT(YEAR FROM cd.date) as year,
        CASE
            WHEN EXTRACT(YEAR FROM cd.date) = fy.first_year THEN 100  -- Use 100 for first year
            ELSE cd.item_index
        END as jan_index
    FROM cpi_data cd
    JOIN first_years fy
        ON cd.item_id = fy.item_id
    WHERE (
        -- For first year, get the first available month
        EXTRACT(YEAR FROM cd.date) = fy.first_year
        AND cd.date = (
            SELECT MIN(date)
            FROM cpi_data cd2
            WHERE EXTRACT(YEAR FROM date) = fy.first_year
                AND cd.item_id = cd2.item_id
        )
        OR EXTRACT(MONTH FROM cd.date) = 1  -- For subsequent years, get January values
    )
),
cumulative_factors AS (
    -- Calculate cumulative product of (index/100) for January values
    SELECT
        year,
        item_id,
        jan_index,
        EXP(
            SUM(LN(jan_index / 100)) OVER (PARTITION BY item_id ORDER BY year)
        ) as cumulative_factor
    FROM january_indices
),
final_data AS (
    -- Join original data with cumulative factors
    SELECT
        cd.date,
        cd.item_id,
        cd.item_index,
        cf.cumulative_factor,
        CASE
            WHEN EXTRACT(MONTH FROM cd.date) = 1
            THEN LAG(cf.cumulative_factor) OVER (PARTITION BY cd.item_id ORDER BY cd.date)
            ELSE cf.cumulative_factor
        END as adjusted_factor
    FROM cpi_data cd
    LEFT JOIN cumulative_factors cf
        ON EXTRACT(YEAR FROM cd.date) = cf.year
        AND cd.item_id = cf.item_id
)

SELECT
    date,
    fd.item_id,
    item_desc,
    item_index AS original_index,
    ROUND(adjusted_factor * item_index, 3) AS rebased_index,
    ROUND(LN(adjusted_factor * item_index), 3) AS ln_rebased_index
FROM final_data fd
    JOIN items i ON fd.item_id = i.item_id
WHERE fd.item_id IN (520206, 520213, 430536, 520249, 520241)
ORDER BY fd.item_id, date;